#!/usr/bin/env python3
"""Fetch gnomAD v4.1 joint frequencies for the MVA genes and the ACMG SF v3.3 genes.

Approved by Santanu in a live chat on 2026-10-04 ("approve gnomAD gene slices"), with the
condition that everything lands inside this project folder and nothing touches C:.

Source: the public gnomAD browser API (https://gnomad.broadinstitute.org/api, dataset
gnomad_r4 = v4.1, joint exome+genome frequencies). Only public gene names and regions are
queried, never the child's variants.
  - 84 ACMG SF v3.3 genes: one gene query each. It returns the gene's exonic variants with a
    small flank, which is all the secondary-findings screen reads.
  - 10 MVA genes: region queries over the whole gene, 60 kb upstream (the known BUB1B
    regulatory second alleles lie ~44 kb upstream) and 5 kb downstream, in chunks of 60 kb or less.
Range queries on the 27-69 GB joint VCFs on S3 were tried first and abandoned: tabix starts each
read at a 16 kb index window, and gnomAD lines are wide, so each exon cost ~10 MB of transfer.

Resource use (the laptop is shared with about 10 daily routines): the process runs at Windows
idle priority. The merge works one chromosome at a time, so peak memory stays small. WSL is
started for bgzip/tabix only if it is already running or at least 1.5 GB of RAM is free.
Otherwise the plain VCF is left in place, the manifest says compression was deferred, and a
rerun (everything is cached by then) finishes the job.

Output (git-ignored): refs/gnomad/gnomad_v4.1_joint_slices.vcf.gz (+ .tbi), sites only, with
INFO AC_joint, AN_joint, AF_joint, nhomalt_joint, fafmax_faf95_max_joint. Raw API responses
are cached in refs/gnomad/api/, so a rerun fetches only what is missing.
Run: python tools/fetch_gnomad_slices.py [--plan-only]
"""
import argparse
import ctypes
import gzip
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from pipeline.sf_screen import SF_V33  # noqa: E402

MANE_GFF = Path(r"D:\mva-hackathon\refs\ensembl\GRCh38.116.mane_select.gff3.gz")
OUT = REPO / "refs" / "gnomad"
API = "https://gnomad.broadinstitute.org/api"
MVA_GENES = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1", "MAD2L1BP", "CEP192", "BUB1"]
ALIASES = {"CENATAC": ["CENATAC", "CCDC84"]}
UPSTREAM, DOWNSTREAM, CHUNK = 60_000, 5_000, 60_000
PAUSE_S = 2.0  # be polite to a free public service
MIN_FREE_RAM_FOR_WSL_GB = 1.5
WSL_TOOLS = "/mnt/d/mva-hackathon/tools/env/bin"
CHROMS = [str(n) for n in range(1, 23)] + ["X", "Y"]
FIELDS = "variant_id pos ref alt joint { ac an homozygote_count fafmax { faf95_max } }"


def lower_priority():
    """Run at Windows idle priority so the shared laptop's other routines keep the CPU."""
    if sys.platform == "win32":
        IDLE_PRIORITY_CLASS = 0x40
        k32 = ctypes.windll.kernel32
        # Declare pointer-sized handles: with ctypes' default int types the pseudo-handle -1 is
        # passed as 0xFFFFFFFF on 64-bit Windows, and the call fails silently.
        k32.GetCurrentProcess.restype = ctypes.c_void_p
        k32.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
        k32.SetPriorityClass.restype = ctypes.c_int
        if not k32.SetPriorityClass(k32.GetCurrentProcess(), IDLE_PRIORITY_CLASS):
            print("warning: could not lower process priority", flush=True)


def free_ram_gb():
    if sys.platform != "win32":
        return None

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

    m = MEMORYSTATUSEX()
    m.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 1e9


def wsl_running():
    try:
        out = subprocess.run(["wsl", "--list", "--running", "--quiet"], capture_output=True, timeout=30).stdout
        return "Ubuntu" in out.decode("utf-16-le", errors="ignore")
    except (OSError, subprocess.TimeoutExpired):
        return False


def mane_genes():
    genes = {}
    with gzip.open(MANE_GFF, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                a = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
                if "Name" in a:
                    genes[a["Name"]] = (f[0], int(f[3]), int(f[4]), f[6])
    return genes


def plan():
    genes = mane_genes()
    jobs, missing = [], []
    for g in MVA_GENES:
        hit = next((n for n in ALIASES.get(g, [g]) if n in genes), None)
        if not hit:
            missing.append(g)
            continue
        chrom, start, end, strand = genes[hit]
        lo, hi = (start - UPSTREAM, end + DOWNSTREAM) if strand == "+" else (start - DOWNSTREAM, end + UPSTREAM)
        lo = max(1, lo)
        for s in range(lo, hi + 1, CHUNK):
            e = min(s + CHUNK - 1, hi)
            jobs.append({"kind": "region", "gene": hit, "chrom": chrom, "start": s, "stop": e,
                         "key": f"region_{hit}_{chrom}_{s}_{e}"})
    for g in sorted(SF_V33):
        if g in genes:
            jobs.append({"kind": "gene", "gene": g, "key": f"gene_{g}"})
        else:
            missing.append(g)
    return jobs, missing


def query(job):
    if job["kind"] == "gene":
        q = f'query {{ gene(gene_symbol: "{job["gene"]}", reference_genome: GRCh38) {{ variants(dataset: gnomad_r4) {{ {FIELDS} }} }} }}'
    else:
        q = (f'query {{ region(chrom: "{job["chrom"]}", start: {job["start"]}, stop: {job["stop"]}, '
             f'reference_genome: GRCh38) {{ variants(dataset: gnomad_r4) {{ {FIELDS} }} }} }}')
    req = urllib.request.Request(API, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": "mva-hackathon-2026"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                body = r.read()
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == 5:
                raise  # a malformed query or unknown gene won't succeed on retry
            time.sleep(min(300, 20 * 2 ** attempt))
            continue
        except (urllib.error.URLError, TimeoutError):
            if attempt == 5:
                raise
            time.sleep(min(300, 20 * 2 ** attempt))
            continue
        d = json.loads(body)
        if d.get("errors"):
            msg = str(d["errors"][0].get("message", d["errors"][0]))
            if "rate" in msg.lower() and attempt < 5:
                time.sleep(30 * (attempt + 1))
                continue
            raise RuntimeError(msg)
        node = d["data"]["gene" if job["kind"] == "gene" else "region"]
        return node["variants"], len(body)
    raise RuntimeError("retries exhausted")


def merge(jobs, cache, plain):
    """Write a sorted, de-duplicated sites VCF, one chromosome at a time."""
    tmp = OUT / "tmp"
    tmp.mkdir(exist_ok=True)
    handles = {}
    try:
        for job in jobs:
            cached = cache / f"{job['key']}.json"
            if not cached.exists():
                continue
            for v in json.loads(cached.read_text(encoding="utf-8")):
                j = v.get("joint")
                if not j or not j.get("an"):
                    continue
                chrom = v["variant_id"].split("-")[0]
                faf = (j.get("fafmax") or {}).get("faf95_max")
                info = (f"AC_joint={j['ac']};AN_joint={j['an']};AF_joint={j['ac'] / j['an']:.4g};"
                        f"nhomalt_joint={j['homozygote_count']}"
                        + (f";fafmax_faf95_max_joint={faf:.4g}" if faf is not None else ""))
                if chrom not in handles:
                    handles[chrom] = open(tmp / f"chr{chrom}.tsv", "w", encoding="ascii", newline="\n")
                handles[chrom].write(f"{v['pos']}\t{v['ref']}\t{v['alt']}\t{info}\n")
    finally:
        for h in handles.values():
            h.close()

    header = [
        "##fileformat=VCFv4.2",
        "##source=gnomAD v4.1 joint (gnomad_r4) via https://gnomad.broadinstitute.org/api",
        '##INFO=<ID=AC_joint,Number=A,Type=Integer,Description="Alternate allele count, exomes+genomes">',
        '##INFO=<ID=AN_joint,Number=1,Type=Integer,Description="Total allele number, exomes+genomes">',
        '##INFO=<ID=AF_joint,Number=A,Type=Float,Description="AC_joint/AN_joint">',
        '##INFO=<ID=nhomalt_joint,Number=A,Type=Integer,Description="Homozygous alternate individuals">',
        '##INFO=<ID=fafmax_faf95_max_joint,Number=A,Type=Float,Description="Max filtering allele frequency (95% CI) across genetic ancestry groups">',
    ] + [f"##contig=<ID=chr{c}>" for c in CHROMS] + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    n = 0
    with open(plain, "w", encoding="ascii", newline="\n") as out:
        out.write("\n".join(header) + "\n")
        for chrom in CHROMS:
            part = tmp / f"chr{chrom}.tsv"
            if not part.exists():
                continue
            rows = part.read_text(encoding="ascii").splitlines()
            rows.sort(key=lambda r: (int(r.split("\t", 1)[0]), r.split("\t")[1], r.split("\t")[2]))
            prev = None
            for r in rows:
                key = r.split("\t", 3)[:3]
                if key == prev:
                    continue  # the same variant returned by two overlapping queries
                prev = key
                pos, ref, alt, info = r.split("\t")
                out.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
                n += 1
            part.unlink()
    tmp.rmdir()
    return n


def wsl_path(p):
    p = str(p)
    return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plan-only", action="store_true", help="print the query plan; fetch nothing")
    args = ap.parse_args()
    lower_priority()
    jobs, missing = plan()
    print(f"{len(jobs)} queries ({sum(j['kind'] == 'region' for j in jobs)} region chunks for MVA genes, "
          f"{sum(j['kind'] == 'gene' for j in jobs)} SF gene queries); missing from MANE: {missing or 'none'}",
          flush=True)
    if args.plan_only:
        return 0

    cache = OUT / "api"
    cache.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    fetched_bytes, failures, per_gene = 0, [], {}
    for i, job in enumerate(jobs, 1):
        cached = cache / f"{job['key']}.json"
        if not cached.exists():
            try:
                variants, n = query(job)
            except Exception as e:  # record and continue; a rerun retries only the missing ones
                failures.append(f"{job['key']}: {type(e).__name__}: {e}")
                print(f"[{i}/{len(jobs)}] {job['key']}: FAILED {e}", flush=True)
                continue
            cached.write_text(json.dumps(variants), encoding="utf-8")
            fetched_bytes += n
            time.sleep(PAUSE_S)
        else:
            variants = json.loads(cached.read_text(encoding="utf-8"))
        per_gene[job["gene"]] = per_gene.get(job["gene"], 0) + len(variants)
        if i % 10 == 0 or i == len(jobs):
            print(f"[{i}/{len(jobs)}] {fetched_bytes / 1e6:.1f} MB fetched this run", flush=True)

    plain = OUT / "gnomad_v4.1_joint_slices.vcf"
    n_rec = merge(jobs, cache, plain)
    print(f"merged {n_rec} records into {plain.name}", flush=True)

    ram = free_ram_gb()
    compressed, note = False, ""
    if wsl_running() or ram is None or ram >= MIN_FREE_RAM_FOR_WSL_GB:
        wp = wsl_path(plain)
        sh = (f"set -e; export HOME=/mnt/d/mva-hackathon/tools/home TMPDIR=/mnt/d/DevCache/tmp/wsl; mkdir -p $TMPDIR; "
              f'nice -n 19 {WSL_TOOLS}/bgzip -f "{wp}" && nice -n 19 {WSL_TOOLS}/tabix -f -p vcf "{wp}.gz"')
        rc = subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", "-c", sh]).returncode
        compressed = rc == 0
        note = "" if compressed else f"bgzip/tabix exit {rc}"
    else:
        note = (f"compression deferred: WSL not running and only {ram:.1f} GB RAM free; "
                "rerun this script (all queries are cached) when the laptop is quiet")
    print(note or "compressed and indexed", flush=True)

    final = OUT / "gnomad_v4.1_joint_slices.vcf.gz"
    manifest = {
        "approved": "Santanu, live chat 2026-10-04: 'approve gnomAD gene slices', inside the project folder, never C:",
        "source": API, "dataset": "gnomad_r4 (v4.1), joint exome+genome",
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "queries": len(jobs), "failures": failures, "missing_from_mane": missing,
        "bytes_fetched_this_run": fetched_bytes, "records": n_rec, "variants_per_gene": per_gene,
        "compressed_and_indexed": compressed, "note": note,
        "output": str(final if compressed else plain),
        "bytes": (final if compressed else plain).stat().st_size if (final if compressed else plain).exists() else None,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("records", "bytes", "bytes_fetched_this_run", "failures",
                                                "compressed_and_indexed")}), flush=True)
    return 1 if failures or (note and not note.startswith("compression deferred")) else 0


if __name__ == "__main__":
    sys.exit(main())
