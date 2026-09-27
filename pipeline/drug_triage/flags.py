"""The six triage flags, as pure functions over a drug record (Open Targets mechanisms + openFDA label text).

Each flag returns {"status", "evidence", "source"} where status is one of:
  "flag"            the source contains a concern for this criterion
  "note"            weaker signal worth a human look (for example, only in the adverse-reactions table)
  "no flag found"   the source was read and nothing matched; this is NOT evidence of safety
  "not assessed"    the source needed for this criterion is missing (for example, no US label)
Paediatric approval uses "yes" / "not established" / "not stated" / "not assessed" instead.

These are text heuristics for triage. Every flag is meant to be read by an expert, not acted on.
"""
import re

# Mitotic-checkpoint and kinetochore genes whose inhibition or degradation weakens the spindle assembly
# checkpoint (Track 2 report §5). For a dosage-deficient query gene, lowering that gene itself is also harmful.
CHECKPOINT_GENES = {"TTK", "AURKB", "INCENP", "PLK1", "BUB1", "BUB1B", "BUB3", "MAD2L1", "MAD1L1",
                    "TRIP13", "CEP57", "CENPE", "ZWINT", "KNL1", "NDC80", "SGO1"}
LOWERING_ACTIONS = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "DEGRADER",
                    "RNAI INHIBITOR", "ANTISENSE INHIBITOR", "DISRUPTING AGENT", "BLOCKER"}

# Mechanism classes that cause chromosome missegregation (aneugens) in their own right.
SPINDLE_GENES = {"AURKA", "KIF11", "KIF18A", "CENPE"}
HDAC_GENES = {f"HDAC{i}" for i in range(1, 12)}

# Split after "." or ";" only when the preceding character is lower-case, a digit, ")" or "]", so that
# abbreviations such as "S. typhimurium" do not cut a negated sentence in two.
_SENT = re.compile(r"(?:(?<=[a-z0-9)\]][.;])|(?<=[a-z0-9)\]] [.;]))\s+(?=[A-Z(0-9])")


def _sentences(text):
    return [s.strip() for s in _SENT.split(text or "") if s.strip()]


def _hit(text, pattern, exclude=None, limit=200):
    """First sentence matching `pattern` (and not `exclude`), trimmed to a window around the match."""
    for s in _sentences(text):
        m = re.search(pattern, s, re.I)
        if m and not (exclude and re.search(exclude, s, re.I)):
            if len(s) <= limit:
                return s
            a = max(0, m.start() - limit // 2)
            return ("..." if a else "") + s[a:a + limit] + ("..." if a + limit < len(s) else "")
    return None


def _res(status, evidence="", source=""):
    return {"status": status, "evidence": evidence, "source": source}


SAFETY_SECTIONS = ("boxed_warning", "warnings_and_cautions", "warnings", "contraindications", "adverse_reactions")


def _unusable(label):
    """A missing label, or one with no safety sections (some unapproved listings), cannot clear a drug."""
    if not label:
        return _res("not assessed", "no US label found", "openFDA")
    if not any(label.get(k) for k in SAFETY_SECTIONS):
        return _res("not assessed", "US label has no warnings, contraindications or adverse-reaction sections", "openFDA")
    return None


def paediatric(label):
    """'yes' if any US prescription label of the ingredient states established paediatric safety and effectiveness."""
    if not label:
        return _res("not assessed", "no US prescription label found", "openFDA")
    texts = label.get("_pediatric_all") or [label.get("pediatric_use", "")]
    texts = [t for t in texts if t]
    if not texts:
        return _res("not assessed", "label has no paediatric-use section", "openFDA")
    for txt in texts:
        s = _hit(txt, r"(safety and (effectiveness|efficacy)|effectiveness and safety)[^.]*\b(ha(ve|s) been|were|was|is) established",
                 exclude=r"\bnot (been )?established|\bnot been\b")
        if s:
            return _res("yes", s, "openFDA 8.4")
    for txt in texts:
        # Only a general statement counts; "not established for migraine" says nothing about other uses.
        s = _hit(txt, r"(safety|effectiveness|efficacy)[^.]*(pediatric|children)[^.]*not (been )?established",
                 exclude=r"not (been )?established (for|in the treatment)")
        if s:
            return _res("not established", s, "openFDA 8.4")
    return _res("not stated", texts[0][:160], "openFDA 8.4")


def checkpoint(drug, query_gene, dosage_deficient=True):
    if not drug or not drug["mechanisms"]:
        return _res("not assessed", "no mechanism of action on record", "Open Targets")
    for m in drug["mechanisms"]:
        lowers = m["action"] in LOWERING_ACTIONS
        hit = set(m["targets"]) & CHECKPOINT_GENES
        if lowers and hit:
            return _res("flag", f"{m['moa']} ({', '.join(sorted(hit))})", "Open Targets MoA")
        if dosage_deficient and lowers and query_gene in m["targets"]:
            return _res("flag", f"lowers {query_gene}, a dosage-deficient protein: {m['moa']}", "Open Targets MoA")
    return _res("no flag found", "; ".join(m["moa"] for m in drug["mechanisms"])[:200], "Open Targets MoA")


def aneugen(drug, label):
    moas = drug["mechanisms"] if drug else []
    for m in moas:
        t = set(m["targets"])
        if any(x.startswith(("TUBA", "TUBB")) for x in t) or re.search(r"tubulin|microtubule", m["moa"], re.I):
            return _res("flag", f"spindle poison class: {m['moa']}", "Open Targets MoA")
        if t & SPINDLE_GENES and m["action"] in LOWERING_ACTIONS:
            return _res("flag", f"mitotic-spindle target: {m['moa']}", "Open Targets MoA")
        if len(t & HDAC_GENES) >= 3 or re.search(r"^histone deacetylase inhibitor$|pan-?HDAC", m["moa"], re.I):
            return _res("flag", f"broad HDAC inhibitor; class causes aneuploidy at high exposure (Olaharski 2006): {m['moa']}",
                        "Open Targets MoA")
    if label:
        tox = " ".join(label.get(k, "") for k in ("nonclinical_toxicology", "carcinogenesis_and_mutagenesis_and_impairment_of_fertility"))
        s = _hit(tox, r"aneugen|aneuploid|polyploid|micronucle|chromosom\w* (aberration|damage|loss)|clastogen",
                 exclude=r"\b(not|negative|no evidence|did not|was not|were not)\b")
        if s:
            return _res("flag", s, "openFDA 13.1")
        if tox:
            return _res("no flag found", "", "openFDA 13.1" + (" + MoA" if moas else ""))
    if moas:
        return _res("no flag found", "", "Open Targets MoA")
    return _res("not assessed", "no mechanism or label", "")


# The drug must be the grammatical subject: "strong CYP3A4 inhibitors increase <drug> exposure" makes it a
# victim, which is handled as a note below.
_PERP = r"(is|was|acts as) an? (strong|moderate) (CYP3A4? )?(inhibitor|inducer)|(inhibits|induces) CYP3A"
_VICTIM = r"substrate of CYP3A|CYP3A4? substrate|metaboli[sz]ed\w* (primarily |mainly |predominantly )?by CYP3A"


def cyp3a4(label):
    if _unusable(label):
        return _unusable(label)
    txt = " ".join(label.get(k, "") for k in ("drug_interactions", "clinical_pharmacology", "pharmacokinetics", "warnings_and_cautions"))
    name = re.escape(label.get("_generic", "").split()[0]) if label.get("_generic") else r"\w+"
    s = _hit(txt, rf"\b{name}\b ({_PERP})", exclude=r"\bnot (an? )?(inhibit|induc)|\bdid not\b|no (clinically )?(significant|relevant)")
    if s:
        return _res("flag", f"perpetrator: {s}", "openFDA 7/12")
    s = _hit(txt, _VICTIM)
    if s:
        return _res("note", f"CYP3A substrate (victim, not perpetrator): {s}", "openFDA 7/12")
    return _res("no flag found", "", "openFDA 7/12")


def kidney(label):
    if _unusable(label):
        return _unusable(label)
    strong = " ".join(label.get(k, "") for k in ("boxed_warning", "contraindications", "warnings_and_cautions", "warnings"))
    s = _hit(strong, r"nephrotox|acute kidney injury|renal (failure|toxicity|tubular)|interstitial nephritis|eGFR (below|less than|<)|contraindicated in patients with (severe )?renal")
    if s:
        return _res("flag", s, "openFDA boxed/W&P/contraindications")
    s = _hit(label.get("adverse_reactions", ""), r"nephrotox|acute kidney injury|renal failure")
    if s:
        return _res("note", s, "openFDA adverse reactions")
    return _res("no flag found", "", "openFDA")


def marrow(label):
    if _unusable(label):
        return _unusable(label)
    pat = r"myelosuppress|bone marrow (suppression|failure|depression)|thrombocytopeni|neutropeni|pancytopeni|agranulocytosis|aplastic an(a)?emia|\bcytopeni"
    ttp = re.compile(r"thrombotic thrombocytopenic purpura", re.I)  # a microangiopathy, not marrow failure
    strong = ttp.sub("TTP", " ".join(label.get(k, "") for k in ("boxed_warning", "warnings_and_cautions", "warnings", "contraindications")))
    s = _hit(strong, pat)
    if s:
        return _res("flag", s, "openFDA boxed/W&P")
    s = _hit(ttp.sub("TTP", label.get("adverse_reactions", "")), pat)
    if s:
        return _res("note", s, "openFDA adverse reactions")
    return _res("no flag found", "", "openFDA")


FLAG_NAMES = ("paediatric_approval", "checkpoint_weakening", "aneugen_class", "cyp3a4_interaction",
              "nephrotoxicity", "myelosuppression")


def evaluate(drug, label, query_gene, dosage_deficient=True):
    return {
        "paediatric_approval": paediatric(label),
        "checkpoint_weakening": checkpoint(drug, query_gene, dosage_deficient),
        "aneugen_class": aneugen(drug, label),
        "cyp3a4_interaction": cyp3a4(label),
        "nephrotoxicity": kidney(label),
        "myelosuppression": marrow(label),
    }


def tier(flags):
    """'excluded' (checkpoint or aneugen), 'caution' (label safety flag), 'screen-clear', or 'insufficient data'."""
    if "flag" in (flags["checkpoint_weakening"]["status"], flags["aneugen_class"]["status"]):
        return "excluded"
    safety = [flags[k]["status"] for k in ("cyp3a4_interaction", "nephrotoxicity", "myelosuppression")]
    if "flag" in safety:
        return "caution"
    if all(s == "not assessed" for s in safety):
        return "insufficient data"
    return "screen-clear"
