#!/usr/bin/env python3
"""Public stability map of BUBR1: predicted effect of every possible amino-acid substitution.

Package R13. Runs as a private Kaggle GPU notebook on public data only: the UniProt sequence and the
AlphaFold DB model of O60566 (BUBR1, 1,050 aa), ESM-2 650M (Lin 2023; MIT) and ThermoMPNN
(Dieckhaus 2024; MIT). No patient data is used; the child's substitution is looked up in the output
locally. How the output is read was fixed beforehand in reports/track2/stability_map_rules.md.

Outputs in /kaggle/working:
  bub1b_stability_map.tsv   pos, wt, mut, esm_llr, ddg_destab (positive = destabilising),
                            ddg_raw, plddt, rsa, contacts
  summary.json              sequence checks, ThermoMPNN sign calibration, control substitutions,
                            the validity gate, and the reference-set percentiles
"""
import glob
import json
import os
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

OUT = Path("/kaggle/working")
WORK = Path("/kaggle/tmp") if Path("/kaggle/tmp").exists() else Path("/tmp")
UNIPROT = "O60566"
AA = "ACDEFGHIKLMNPQRSTVWY"
# Suijkerbuijk 2010 (PMID 20516114): tested substitutions and their class.
CONTROLS = {"R727C": "abundance", "R814H": "abundance", "L844F": "abundance", "I909T": "abundance",
            "L1012P": "abundance", "R550Q": "no defect", "Q921H": "no defect", "Y155C": "functional"}
# Theoretical maximum ASA (Tien 2013), for relative solvent accessibility.
MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, E=223, Q=225, G=104, H=224, I=197, L=201, K=236,
               M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
T0 = time.time()
summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def fetch(url, binary=False, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read()
                return data if binary else data.decode()
        except Exception as e:  # network hiccup: back off and retry
            log("fetch failed", i + 1, url, e)
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"could not fetch {url}")


def sh(cmd, check=True):
    log("$", cmd)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    tail = (p.stdout + p.stderr)[-3000:]
    if tail.strip():
        print(tail, flush=True)
    if check and p.returncode:
        raise RuntimeError(f"failed ({p.returncode}): {cmd}")
    return p


def save():
    summary["seconds"] = round(time.time() - T0, 1)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))


# ---------------------------------------------------------------- sequence and structure
seq = "".join(fetch(f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta").split("\n")[1:]).strip()
L = len(seq)
log("UniProt length", L)
pred = json.loads(fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}"))[0]
pdb_path = WORK / f"AF-{UNIPROT}.pdb"
pdb_path.write_bytes(fetch(pred["pdbUrl"], binary=True))
summary["alphafold"] = {"url": pred["pdbUrl"], "version": pred.get("latestVersion"), "sequence_matches_uniprot": pred.get("uniprotSequence") == seq}
summary["controls_wt_ok"] = {k: seq[int(k[1:-1]) - 1] == k[0] for k in CONTROLS}
log("controls wt ok", summary["controls_wt_ok"])
assert all(summary["controls_wt_ok"].values()), "control residues do not match the UniProt sequence"

sh(f"{sys.executable} -m pip install -q biopython fair-esm", check=False)
from Bio.PDB import PDBParser  # noqa: E402
from Bio.PDB.SASA import ShrakeRupley  # noqa: E402

structure = PDBParser(QUIET=True).get_structure("bubr1", str(pdb_path))
chain = next(structure[0].get_chains())
ShrakeRupley().compute(structure[0], level="R")
residues = {r.id[1]: r for r in chain if r.id[0] == " "}
plddt = np.full(L + 1, np.nan)
rsa = np.full(L + 1, np.nan)
cb = {}
for i, r in residues.items():
    plddt[i] = next(iter(r)).get_bfactor()
    rsa[i] = min(1.0, r.sasa / MAX_ASA[seq[i - 1]])
    atom = r["CB"] if "CB" in r else r["CA"]
    cb[i] = atom.coord
coords = np.array([cb[i] for i in range(1, L + 1)])
dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
contacts = np.concatenate([[np.nan], (dist < 10).sum(1) - 1])
log("structure features done; mean pLDDT", float(np.nanmean(plddt[1:])))

# ---------------------------------------------------------------- ESM-2 masked marginals
import torch  # noqa: E402
import esm  # noqa: E402

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
dev = "cuda" if torch.cuda.is_available() else "cpu"
model = model.eval().to(dev)
if dev == "cuda":
    model = model.half()
summary["esm"] = {"model": "esm2_t33_650M_UR50D", "device": dev}
conv = alphabet.get_batch_converter()
W = 1022
aa_idx = [alphabet.get_idx(a) for a in AA]
llr = np.full((L + 1, 20), np.nan)
B = 8 if dev == "cuda" else 1
with torch.no_grad():
    for start in range(1, L + 1, B):
        positions = list(range(start, min(L, start + B - 1) + 1))
        batch, offs = [], []
        for p in positions:
            s = min(max(0, p - 1 - W // 2), max(0, L - W))
            batch.append((str(p), seq[s:s + W]))
            offs.append(p - 1 - s + 1)  # +1 for BOS
        _, _, toks = conv(batch)
        for k, o in enumerate(offs):
            toks[k, o] = alphabet.mask_idx
        logits = model(toks.to(dev))["logits"].float()
        lp = torch.log_softmax(logits, -1)
        for k, (p, o) in enumerate(zip(positions, offs)):
            row = lp[k, o, aa_idx].cpu().numpy()
            llr[p] = row - row[AA.index(seq[p - 1])]
        if start % 200 == 1:
            log("ESM position", start)
summary["esm"]["positions"] = int(np.isfinite(llr[1:, 0]).sum())
save()

# ---------------------------------------------------------------- ThermoMPNN
ddg_raw = np.full((L + 1, 20), np.nan)
tm = WORK / "ThermoMPNN"
try:
    sh(f"git clone -q --depth 1 https://github.com/Kuhlman-Lab/ThermoMPNN {tm}")
    sh(f"{sys.executable} -m pip install -q omegaconf pytorch-lightning", check=False)
    commit = sh(f"git -C {tm} rev-parse HEAD").stdout.strip()
    summary["thermompnn"] = {"commit": commit}
    sh(f"ls -R {tm} | head -60", check=False)
    script = next(iter(glob.glob(f"{tm}/analysis/custom_inference.py") + glob.glob(f"{tm}/**/custom_inference.py", recursive=True)), None)
    if not script:
        raise RuntimeError("custom_inference.py not found")
    sh(f"sed -n '1,200p' {script} | grep -n 'add_argument\\|yaml\\|load(' ", check=False)
    local_yaml = Path(tm) / "local.yaml"
    if not local_yaml.exists():
        local_yaml.write_text(f"platform:\n  accel: {'gpu' if dev == 'cuda' else 'cpu'}\n  thermompnn_dir: {tm}\n")
    outdir = WORK / "tm_out"
    outdir.mkdir(exist_ok=True)
    model_path = next(iter(glob.glob(f"{tm}/models/thermoMPNN_default.pt") + glob.glob(f"{tm}/**/*default*.pt", recursive=True)))
    run = sh(f"cd {Path(script).parent} && {sys.executable} {script} --pdb {pdb_path} --chain A "
             f"--model_path {model_path} --out_dir {outdir}", check=False)
    csvs = glob.glob(f"{outdir}/**/*.csv", recursive=True) + glob.glob(f"{Path(script).parent}/*.csv")
    if not csvs:
        raise RuntimeError("ThermoMPNN wrote no CSV")
    import pandas as pd  # noqa: E402
    df = pd.read_csv(max(csvs, key=os.path.getmtime))
    log("ThermoMPNN columns", list(df.columns), len(df))
    mcol = next(c for c in df.columns if df[c].astype(str).str.match(r"^[A-Z]\d+[A-Z]$").mean() > 0.5)
    dcol = next(c for c in df.columns if "ddg" in c.lower())
    for m, v in zip(df[mcol].astype(str), df[dcol]):
        mm = re.match(r"^([A-Z])(\d+)([A-Z])$", m)
        if mm and mm.group(3) in AA:
            p = int(mm.group(2))
            if 1 <= p <= L and seq[p - 1] == mm.group(1):
                ddg_raw[p, AA.index(mm.group(3))] = float(v)
    summary["thermompnn"].update({"file": os.path.basename(max(csvs, key=os.path.getmtime)), "mutation_col": mcol,
                                  "ddg_col": dcol, "scored": int(np.isfinite(ddg_raw[1:]).sum())})
except Exception as e:  # keep the ESM and structure results even if ThermoMPNN fails
    summary["thermompnn_error"] = repr(e)
    log("ThermoMPNN failed:", e)
save()

# ---------------------------------------------------------------- sign calibration
sign = 1.0
if np.isfinite(ddg_raw).any():
    buried = [p for p in range(1, L + 1) if rsa[p] < 0.1 and plddt[p] >= 70 and seq[p - 1] in "LIVFMW"]
    exposed = [p for p in range(1, L + 1) if rsa[p] > 0.5 and plddt[p] >= 70]
    harsh = [AA.index(a) for a in "DEKRP"]
    m_buried = float(np.nanmean(ddg_raw[np.ix_(buried, harsh)]))
    m_exposed = float(np.nanmean(ddg_raw[exposed, AA.index("A")]))
    sign = 1.0 if m_buried > m_exposed else -1.0
    summary["sign_calibration"] = {"buried_hydrophobic_to_DEKRP_mean": m_buried, "exposed_to_A_mean": m_exposed,
                                   "n_buried": len(buried), "n_exposed": len(exposed),
                                   "destabilising_is": "positive raw" if sign > 0 else "negative raw"}
ddg = sign * ddg_raw

# ---------------------------------------------------------------- table, controls, gate
with open(OUT / "bub1b_stability_map.tsv", "w") as fh:
    fh.write("pos\twt\tmut\tesm_llr\tddg_destab\tddg_raw\tplddt\trsa\tcontacts\n")
    for p in range(1, L + 1):
        for j, a in enumerate(AA):
            if a == seq[p - 1]:
                continue
            fh.write(f"{p}\t{seq[p - 1]}\t{a}\t{llr[p, j]:.4f}\t{ddg[p, j]:.4f}\t{ddg_raw[p, j]:.4f}\t"
                     f"{plddt[p]:.1f}\t{rsa[p]:.3f}\t{contacts[p]:.0f}\n")

ref = [(p, j) for p in range(727, L + 1) if plddt[p] >= 70 for j, a in enumerate(AA) if a != seq[p - 1]]
ref_ddg = np.array([ddg[p, j] for p, j in ref])
ref_llr = np.array([llr[p, j] for p, j in ref])


def pct(values, v):
    values = values[np.isfinite(values)]
    return float((values < v).mean() * 100) if len(values) and np.isfinite(v) else None


ctl = {}
for name, cls in CONTROLS.items():
    p, j = int(name[1:-1]), AA.index(name[-1])
    ctl[name] = {"class": cls, "ddg_destab": float(ddg[p, j]), "esm_llr": float(llr[p, j]), "plddt": float(plddt[p]),
                 "rsa": float(rsa[p]), "ddg_pct_ref": pct(ref_ddg, ddg[p, j]), "esm_pct_ref": pct(-ref_llr, -llr[p, j])}
summary["controls"] = ctl
summary["reference_set"] = {"definition": "all substitutions at residues 727-1050 with pLDDT >= 70", "n": len(ref),
                            "ddg_p25": float(np.nanpercentile(ref_ddg, 25)) if np.isfinite(ref_ddg).any() else None,
                            "ddg_p75": float(np.nanpercentile(ref_ddg, 75)) if np.isfinite(ref_ddg).any() else None}
if np.isfinite(ref_ddg).any():
    nd = max(ctl["R550Q"]["ddg_destab"], ctl["Q921H"]["ddg_destab"])
    above = sum(ctl[k]["ddg_destab"] > nd for k in ("R727C", "R814H", "L844F", "I909T", "L1012P"))
    summary["gate"] = {"abundance_controls_above_both_no_defect": above, "passes": above >= 3}
nd_llr = min(ctl["R550Q"]["esm_llr"], ctl["Q921H"]["esm_llr"])
summary["esm_controls_below_both_no_defect"] = sum(ctl[k]["esm_llr"] < nd_llr for k in ("R727C", "R814H", "L844F", "I909T", "L1012P"))
summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
save()
log(json.dumps({k: summary.get(k) for k in ("gate", "sign_calibration", "thermompnn_error")}, indent=1))
log("done")
