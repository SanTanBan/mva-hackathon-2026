#!/usr/bin/env python3
"""Public splice map of BUB1B: Pangolin scores for every possible single-nucleotide change.

Package R09. Runs as a private Kaggle GPU notebook on public data only: the GRCh38 sequence of the
BUB1B locus and the RefSeq transcript NM_001211.6 (MANE Select), both fetched from NCBI, and the
Pangolin model (Zeng & Li 2022, Genome Biol 23:103; GPL-3.0), cloned from GitHub. No patient data is
used. The child's variants are later looked up in this table locally, so they never leave the
approved environment.

Scoring follows Pangolin's own compute_score: for each tissue, average (alt - ref) over its models,
take the per-position minimum (loss) and maximum (gain) across tissues within +/-50 nt, and mask
annotated splice sites (gain at annotated sites set to <= 0, loss at unannotated sites set to >= 0).
Masking uses MANE Select exon boundaries only. The batched implementation is checked against
Pangolin's compute_score on a sample of variants before the full run.
"""
import json
import os
import re
import subprocess
import sys
import threading
import time
import types
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import torch

OUT = Path("/kaggle/working")
WORK = Path("/kaggle/tmp") if Path("/kaggle/tmp").exists() else Path("/tmp")
CHROM_ACC, CHROM = "NC_000015.10", "chr15"
TX = "NM_001211.6"
REGION = (40_150_000, 40_232_000)   # 1-based inclusive; covers BUB1B with >= 5,050 nt of flank
D = 50                               # Pangolin's default scoring distance
FLANK = 5000
BATCH = int(os.environ.get("BATCH", "192"))  # variants per batch (64 positions x 3)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def fetch(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # network hiccup: back off and retry
            log("fetch failed", i + 1, e)
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"could not fetch {url}")


# ---------------------------------------------------------------- reference and transcript
fa = fetch(f"{EUTILS}?db=nuccore&id={CHROM_ACC}&seq_start={REGION[0]}&seq_stop={REGION[1]}&rettype=fasta&retmode=text")
region_seq = "".join(fa.split("\n")[1:]).upper()
assert len(region_seq) == REGION[1] - REGION[0] + 1, len(region_seq)
gb = fetch(f"{EUTILS}?db=nuccore&id={TX}&rettype=gb&retmode=text")
mrna = "".join(re.findall(r"[acgtn]+", gb.split("ORIGIN", 1)[1])).upper()
exons_m = [tuple(map(int, m)) for m in re.findall(r"^     exon\s+(\d+)\.\.(\d+)", gb, re.M)]
cds_s, cds_e = map(int, re.search(r"^     CDS\s+(\d+)\.\.(\d+)", gb, re.M).groups())
log("region", len(region_seq), "nt; mRNA", len(mrna), "nt;", len(exons_m), "exons; CDS", cds_s, cds_e)


def gpos(i):
    return REGION[0] + i


exons = []  # (exon number, genomic start, genomic end, mRNA start, mRNA end)
for k, (ms, me) in enumerate(exons_m, 1):
    seg = mrna[ms - 1:me]
    idx = region_seq.find(seg)
    assert idx >= 0 and region_seq.find(seg, idx + 1) < 0, f"exon {k} not uniquely placed"
    exons.append((k, gpos(idx), gpos(idx) + len(seg) - 1, ms, me))
assert all(a[2] < b[1] for a, b in zip(exons, exons[1:])), "exons out of order (strand?)"
TX_START, TX_END = exons[0][1], exons[-1][2]
log("MANE span", CHROM, TX_START, TX_END, "=", TX_END - TX_START + 1, "nt")
ANNOTATED = sorted({e[1] for e in exons} | {e[2] for e in exons})


def c_of_mrna(m):
    if m < cds_s:
        return f"-{cds_s - m}"
    if m > cds_e:
        return f"*{m - cds_e}"
    return str(m - cds_s + 1)


def hgvs_c(p):
    """HGVS c. position of a genomic position within the MANE transcript (plus strand)."""
    for k, gs, ge, ms, me in exons:
        if gs <= p <= ge:
            return c_of_mrna(ms + p - gs), "exon", k
    for (k, gs, ge, ms, me), (k2, gs2, ge2, ms2, me2) in zip(exons, exons[1:]):
        if ge < p < gs2:
            up, down = p - ge, gs2 - p
            if up <= down:
                return f"{c_of_mrna(me)}+{up}", "intron", k
            return f"{c_of_mrna(ms2)}-{down}", "intron", k
    raise ValueError(p)


# ---------------------------------------------------------------- Pangolin
os.chdir(WORK)  # clone outside /kaggle/working so Pangolin itself is not saved as notebook output
if not Path("Pangolin").exists():
    subprocess.run(["git", "clone", "--depth", "1", "https://github.com/tkzeng/Pangolin.git"], check=True)
commit = subprocess.run(["git", "-C", "Pangolin", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
sys.path.insert(0, str(Path("Pangolin").resolve()))
import pangolin.model as pm  # noqa: E402

files = sorted(os.listdir("Pangolin/pangolin/models"))
v2 = any(f.endswith(".v2") for f in files)
pat = re.compile(r"^final\.(\d+)\.(\d+)\.3" + (r"\.v2$" if v2 else r"$"))
reps = sorted({int(m.group(1)) for f in files if (m := pat.match(f))})
src = Path("Pangolin/pangolin/pangolin.py").read_text()
per_tissue = int(re.search(r"models\[(\d+)\*j:", src).group(1))
log("Pangolin commit", commit[:10], "| v2 weights:", v2, "| replicates", reps, "| models per tissue in source:", per_tissue)
if per_tissue != len(reps):
    reps = reps[:per_tissue]
TISSUE_CH = [1, 4, 7, 10]


def load_models(device):
    models = []
    for i in [0, 2, 4, 6]:
        for j in reps:
            m = pm.Pangolin(pm.L, pm.W, pm.AR)
            name = f"final.{j}.{i}.3" + (".v2" if v2 else "")
            w = torch.load(f"Pangolin/pangolin/models/{name}", map_location=device)
            m.load_state_dict(w)
            m.to(device).eval()
            models.append(m)
    return models


torch.backends.cudnn.benchmark = True
NGPU = torch.cuda.device_count()
devices = [torch.device(f"cuda:{i}") for i in range(NGPU)] or [torch.device("cpu")]
log("devices", devices, [torch.cuda.get_device_name(i) for i in range(NGPU)])
MODELS = {d: load_models(d) for d in devices}
N_PER = len(reps)

LUT = np.zeros(256, dtype=np.int64)
for ch, v in zip("ACGT", range(1, 5)):
    LUT[ord(ch)] = v
EYE = np.vstack([np.zeros(4, dtype=np.float32), np.eye(4, dtype=np.float32)])


def onehot(seqs):
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)
    return EYE[LUT[arr]].transpose(0, 2, 1)  # (B, 4, len)


@torch.no_grad()
def tissue_preds(x, device, amp):
    """x: (B, 4, len) -> (4 tissues, B, len - 10000): mean P(splice) over each tissue's models."""
    xt = torch.from_numpy(np.ascontiguousarray(x)).to(device)
    out = []
    with torch.autocast("cuda", dtype=torch.float16, enabled=amp and device.type == "cuda"):
        for t in range(4):
            acc = None
            for m in MODELS[device][t * N_PER:(t + 1) * N_PER]:
                y = m(xt)[:, TISSUE_CH[t], :].float()
                acc = y if acc is None else acc + y
            out.append((acc / N_PER).cpu().numpy())
    return np.stack(out)


scan = list(range(TX_START, TX_END + 1))
ann_lo, ann_hi = TX_START - D, TX_END + D
ANN = np.zeros(ann_hi - ann_lo + 1, dtype=bool)
for a in ANNOTATED:
    if ann_lo <= a <= ann_hi:
        ANN[a - ann_lo] = True


def ref_window(p):
    i = p - REGION[0]
    return region_seq[i - FLANK - D:i + FLANK + D + 1]


def score_variants(variants, device, amp=True):
    """variants: list of (pos, alt). Each variant is scored against a reference window computed in the
    same batch and precision, so rounding errors largely cancel. Returns Pangolin-style gain/loss rows."""
    ref_idx, windows = {}, []
    for p, _ in variants:
        if p not in ref_idx:
            ref_idx[p] = len(windows)
            windows.append(ref_window(p))
    alt_idx = []
    for p, a in variants:
        r = windows[ref_idx[p]]
        alt_idx.append(len(windows))
        windows.append(r[:FLANK + D] + a + r[FLANK + D + 1:])
    pred = tissue_preds(onehot(windows), device, amp)   # (4, n_windows, 2D+1)
    rows = []
    for (p, a), ai in zip(variants, alt_idx):
        diff = pred[:, ai, :] - pred[:, ref_idx[p], :]
        loss = diff[np.argmin(diff, axis=0), np.arange(2 * D + 1)]
        gain = diff[np.argmax(diff, axis=0), np.arange(2 * D + 1)]
        k = p - D - ann_lo
        ann = ANN[k:k + 2 * D + 1]
        gain = np.where(ann, np.minimum(gain, 0), gain)
        loss = np.where(~ann, np.maximum(loss, 0), loss)
        g, l = int(np.argmax(gain)), int(np.argmin(loss))
        rows.append((p, a, round(float(gain[g]), 4), g - D, round(float(loss[l]), 4), l - D))
    return rows


# ---------------------------------------------------------------- check against Pangolin's own scorer
for mod in ("vcf", "gffutils", "pyfastx"):
    sys.modules.setdefault(mod, types.ModuleType(mod))
try:
    import pkg_resources  # noqa: F401  (pangolin.py imports it at module level)
except ImportError:
    _pr = types.ModuleType("pkg_resources")
    _pr.resource_filename = lambda *a: ""
    sys.modules["pkg_resources"] = _pr
check = {"compared": 0, "max_abs_diff_fp32": None, "max_abs_diff_amp": None, "note": ""}
try:
    import pangolin.pangolin as pp
    rng = np.random.default_rng(7)
    sample = [(int(p), "ACGT"[int(rng.integers(4))]) for p in rng.choice(scan, 12, replace=False)]
    sample = [(p, a if a != region_seq[p - REGION[0]] else "T" if region_seq[p - REGION[0]] != "T" else "A") for p, a in sample]
    mine32 = {(r[0], r[1]): r for r in score_variants(sample, devices[0], amp=False)}
    mine16 = {(r[0], r[1]): r for r in score_variants(sample, devices[0], amp=True)}
    diffs, diffs16 = [], []
    for p, a in sample:
        i = p - REGION[0]
        ref_s = region_seq[i - FLANK - D:i + FLANK + D + 1]
        alt_s = ref_s[:FLANK + D] + a + ref_s[FLANK + D + 1:]
        loss, gain = pp.compute_score(ref_s, alt_s, "+", D, MODELS[devices[0]])
        pos = np.array([x for x in ANNOTATED]) - (p - D)
        pos = pos[(pos >= 0) & (pos < len(loss))]
        gain[pos] = np.minimum(gain[pos], 0)
        notp = ~np.isin(np.arange(len(loss)), pos)
        loss[notp] = np.maximum(loss[notp], 0)
        diffs += [abs(float(np.max(gain)) - mine32[(p, a)][2]), abs(float(np.min(loss)) - mine32[(p, a)][4])]
        diffs16 += [abs(float(np.max(gain)) - mine16[(p, a)][2]), abs(float(np.min(loss)) - mine16[(p, a)][4])]
    check.update(compared=len(sample), max_abs_diff_fp32=max(diffs), max_abs_diff_amp=max(diffs16))
    log("check vs Pangolin compute_score:", check)
    assert max(diffs) < 1e-3, "batched scorer disagrees with Pangolin"
    assert max(diffs16) < 2e-2, "mixed precision drifts too far from Pangolin"
except AssertionError:
    raise
except Exception as e:  # scorer import failed: record it and continue with the batched scorer
    check["note"] = f"could not import Pangolin's scorer: {type(e).__name__}: {e}"
    log(check["note"])

# ---------------------------------------------------------------- positive controls (published ClinVar/literature variants)
CONTROLS = {"c.1402-1G>T": "loss", "c.2386-2A>G": "loss", "c.581+1G>T": "loss", "c.1401G>T": "loss",
            "c.1288+5G>A": "loss", "c.1402-5A>G": "either", "c.2386-11A>G": "gain"}
cmap = {}
for p in scan:
    cmap[hgvs_c(p)[0]] = p
ctrl_vars = []
for h, _ in CONTROLS.items():
    m = re.match(r"c\.([-*]?\d+(?:[+-]\d+)?)([ACGT])>([ACGT])", h)
    p = cmap[m.group(1)]
    assert region_seq[p - REGION[0]] == m.group(2), f"{h}: reference base mismatch, mapping is wrong"
    ctrl_vars.append((h, p, m.group(3)))
ctrl_rows = score_variants([(p, a) for _, p, a in ctrl_vars], devices[0], amp=True)
ctrl_rows32 = score_variants([(p, a) for _, p, a in ctrl_vars], devices[0], amp=False)
controls = []
for (h, p, a), r, r32 in zip(ctrl_vars, ctrl_rows, ctrl_rows32):
    controls.append({"hgvs_c": h, "pos": p, "gain": r[2], "gain_at": r[3], "loss": r[4], "loss_at": r[5],
                     "gain_fp32": r32[2], "loss_fp32": r32[4], "expected": CONTROLS[h]})
    log("control", h, "gain", r[2], "@", r[3], "loss", r[4], "@", r[5])

# ---------------------------------------------------------------- full saturation scan, split across GPUs
variants = [(p, a) for p in scan for a in "ACGT" if a != region_seq[p - REGION[0]]]
log("variants to score:", len(variants), "| batch", BATCH, "| GPUs", len(devices))
# shard by position so the three alternatives at a position share one reference window
shards = [[v for v in variants if (v[0] - TX_START) % len(devices) == k] for k in range(len(devices))]
results = [[] for _ in devices]
progress = [0] * len(devices)


def run_shard(k):
    dev, vs = devices[k], shards[k]
    for i in range(0, len(vs), BATCH):
        results[k].extend(score_variants(vs[i:i + BATCH], dev, amp=True))
        progress[k] = min(i + BATCH, len(vs))
        if (i // BATCH) % 100 == 0:
            log(f"gpu{k}", progress[k], "/", len(vs))


with ThreadPoolExecutor(len(devices)) as ex:
    list(ex.map(run_shard, range(len(devices))))
rows = sorted(r for rs in results for r in rs)
log("scored", len(rows))

# ---------------------------------------------------------------- write outputs
import gzip  # noqa: E402

with gzip.open(OUT / "bub1b_pangolin_snv.tsv.gz", "wt") as fh:
    fh.write("chrom\tpos\tref\talt\thgvs_c\tregion\texon_or_intron\tgain\tgain_offset\tloss\tloss_offset\n")
    for p, a, g, go, l, lo in rows:
        h, reg, num = hgvs_c(p)
        fh.write(f"{CHROM}\t{p}\t{region_seq[p - REGION[0]]}\t{a}\tc.{h}{region_seq[p - REGION[0]]}>{a}\t{reg}\t{num}\t{g}\t{go}\t{l}\t{lo}\n")

meta = {
    "transcript": TX, "chrom_accession": CHROM_ACC, "region": REGION, "mane_span": [TX_START, TX_END],
    "exons": [dict(zip(["exon", "g_start", "g_end", "m_start", "m_end"], e)) for e in exons],
    "pangolin_commit": commit, "pangolin_v2_weights": v2, "models_per_tissue": N_PER, "distance": D,
    "masking": "MANE Select exon boundaries only", "precision": "float16 autocast; each variant scored against a reference window from the same batch", "n_variants": len(rows), "check_vs_pangolin": check,
    "controls": controls, "gpus": [torch.cuda.get_device_name(i) for i in range(NGPU)],
    "runtime_s": round(time.time() - T0), "torch": torch.__version__,
}
(OUT / "bub1b_pangolin_meta.json").write_text(json.dumps(meta, indent=1))
log("done")
