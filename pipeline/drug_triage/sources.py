"""Public-data fetchers for the drug-triage pipeline: Open Targets, ChEMBL and openFDA drug labels.

Only public identifiers leave this machine (gene symbols, drug names, ChEMBL and Ensembl IDs).
Nothing patient-derived is ever passed to these functions. Responses are cached as JSON outside the
repository so reruns are reproducible and polite to the APIs.
"""
import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"
FDA_URL = "https://api.fda.gov/drug/label.json"
UA = "mva-hackathon-drug-triage/0.1 (public-data research)"

DEFAULT_CACHE = Path(os.environ.get("DRUG_TRIAGE_CACHE", r"D:\mva-hackathon\refs\drug_triage_cache"))

# openFDA label sections the flags read. Everything else is dropped before caching.
LABEL_SECTIONS = (
    "boxed_warning", "contraindications", "warnings_and_cautions", "warnings", "precautions",
    "adverse_reactions", "drug_interactions", "clinical_pharmacology", "pharmacokinetics",
    "pediatric_use", "use_in_specific_populations", "nonclinical_toxicology",
    "carcinogenesis_and_mutagenesis_and_impairment_of_fertility", "indications_and_usage",
)


class Http:
    """Cached JSON over HTTP. `offline=True` raises on a cache miss (used by tests and dry runs)."""

    def __init__(self, cache_dir=DEFAULT_CACHE, offline=False, pause=0.25):
        self.cache_dir = Path(cache_dir)
        self.offline = offline
        self.pause = pause

    def _key(self, method, url, body):
        h = hashlib.sha256(f"{method} {url} {json.dumps(body, sort_keys=True)}".encode()).hexdigest()
        return self.cache_dir / f"{h[:2]}" / f"{h}.json"

    def request(self, url, body=None):
        method = "POST" if body is not None else "GET"
        path = self._key(method, url, body)
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
        if self.offline:
            raise LookupError(f"cache miss in offline mode: {url}")
        data = json.dumps(body).encode() if body is not None else None
        headers = {"Accept": "application/json", "User-Agent": UA}
        if data:
            headers["Content-Type"] = "application/json"
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    payload = json.load(r)
                break
            except urllib.error.HTTPError as e:
                if e.code == 404:  # openFDA answers "no match" with 404
                    payload = {"_not_found": True}
                    break
                if e.code in (429, 500, 502, 503) and attempt < 2:
                    time.sleep(2 * (attempt + 1))
                    continue
                raise
        time.sleep(self.pause)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return payload


def _ot(http, query, variables):
    out = http.request(OT_URL, {"query": query, "variables": variables})
    if out.get("errors"):
        raise RuntimeError(f"Open Targets error: {out['errors'][0].get('message')}")
    return out["data"]


def ot_target(http, symbol):
    """Ensembl ID, tractability flags and drug/clinical candidates for a gene symbol."""
    hits = _ot(http, "query($q:String!){ search(queryString:$q, entityNames:[\"target\"]){ hits { id name } } }",
               {"q": symbol})["search"]["hits"]
    ens = next((h["id"] for h in hits if h["name"].upper() == symbol.upper()), None)
    if not ens:
        return {"symbol": symbol, "ensembl": None, "tractability": [], "candidates": []}
    t = _ot(http, """query($id:String!){ target(ensemblId:$id){ approvedSymbol
            tractability { label modality value }
            drugAndClinicalCandidates { count rows { maxClinicalStage drug { id name } } } } }""", {"id": ens})["target"]
    rows = (t.get("drugAndClinicalCandidates") or {}).get("rows") or []
    return {
        "symbol": t["approvedSymbol"], "ensembl": ens,
        "tractability": [x for x in t.get("tractability") or [] if x["value"]],
        "candidates": [{"chembl": r["drug"]["id"], "name": r["drug"]["name"], "stage": r["maxClinicalStage"]} for r in rows],
    }


def ot_drug(http, name=None, chembl=None):
    """Mechanisms, maximum clinical stage and warnings for a drug, looked up by ChEMBL ID or name."""
    if not chembl:
        hits = _ot(http, "query($q:String!){ search(queryString:$q, entityNames:[\"drug\"]){ hits { id name } } }",
                   {"q": name})["search"]["hits"]
        chembl = next((h["id"] for h in hits if h["name"].upper() == name.upper()), None) or (hits[0]["id"] if hits else None)
    if not chembl:
        return None
    d = _ot(http, """query($id:String!){ drug(chemblId:$id){ id name maximumClinicalStage drugType
            mechanismsOfAction { rows { mechanismOfAction actionType targets { approvedSymbol } } }
            drugWarnings { warningType toxicityClass } parentMolecule { id } } }""", {"id": chembl})["drug"]
    if not d:
        return None
    return {
        "chembl": d["id"], "parent": (d.get("parentMolecule") or {}).get("id") or d["id"], "name": d["name"], "stage": d.get("maximumClinicalStage"), "type": d.get("drugType"),
        "mechanisms": [{"moa": r["mechanismOfAction"], "action": (r.get("actionType") or "").upper(),
                        "targets": [t["approvedSymbol"] for t in r.get("targets") or []]}
                       for r in (d.get("mechanismsOfAction") or {}).get("rows") or []],
        "warnings": [{"type": w.get("warningType"), "class": w.get("toxicityClass")} for w in d.get("drugWarnings") or []],
    }


def chembl_potency(http, symbol):
    """Best pChEMBL and count of pChEMBL >= 6 records against the human single-protein target for a symbol."""
    res = http.request(f"{CHEMBL_URL}/target/search.json?" + urllib.parse.urlencode({"q": symbol, "limit": 20}))
    tid = None
    for t in res.get("targets", []):
        if t.get("organism") != "Homo sapiens" or t.get("target_type") != "SINGLE PROTEIN":
            continue
        syns = {s["component_synonym"].upper() for c in t.get("target_components", []) for s in c.get("target_component_synonyms", [])}
        if symbol.upper() in syns:
            tid = t["target_chembl_id"]
            break
    if not tid:
        return {"target": None, "records": 0, "potent": 0, "best_pchembl": None}
    act = http.request(f"{CHEMBL_URL}/activity.json?" + urllib.parse.urlencode(
        {"target_chembl_id": tid, "pchembl_value__isnull": "false", "only": "pchembl_value", "limit": 1000}))
    vals = [float(a["pchembl_value"]) for a in act.get("activities", []) if a.get("pchembl_value")]
    return {"target": tid, "records": act.get("page_meta", {}).get("total_count", len(vals)),
            "potent": sum(v >= 6 for v in vals), "best_pchembl": max(vals) if vals else None}


def fda_label(http, name):
    """The most recent US prescription label whose generic name matches, reduced to the sections the flags read.

    OTC and homeopathic listings are skipped. `_pediatric_all` keeps the paediatric-use text of every matching
    label, because generics and formulations of one active ingredient often differ in what they state.
    """
    q = urllib.parse.urlencode({"search": f'openfda.generic_name:"{name.lower()}" AND openfda.product_type:"HUMAN PRESCRIPTION DRUG"',
                                "limit": 20})
    res = http.request(f"{FDA_URL}?{q}")
    # Single-ingredient labels only: "METFORMIN HYDROCHLORIDE" matches metformin, while combination products and
    # multi-ingredient homeopathic listings that merely contain the name do not.
    want = name.upper()
    results = [r for r in res.get("results") or []
               if any(g.upper() == want or (g.upper().startswith(want + " ") and " AND " not in g.upper() and "," not in g)
                      for g in (r.get("openfda") or {}).get("generic_name", []))]
    if not results:
        return None
    lab = max(results, key=lambda r: r.get("effective_time", ""))
    out = {k: " ".join(lab[k]) for k in LABEL_SECTIONS if k in lab}
    out["_set_id"] = lab.get("set_id")
    out["_effective"] = lab.get("effective_time")
    out["_generic"] = (lab.get("openfda") or {}).get("generic_name", [name])[0]
    out["_pediatric_all"] = [" ".join(r["pediatric_use"]) for r in results if r.get("pediatric_use")]
    out["_n_labels"] = len(results)
    return out
