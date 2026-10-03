#!/usr/bin/env python3
"""Fetch gnomAD v4.1 joint-sites slices for the MVA genes and the ACMG SF v3.3 genes.

Approved by Santanu in a live chat on 2026-10-04 ("approve gnomAD gene slices"), with the
condition that everything lands inside this project folder and nothing touches C:.

Only public gene regions are queried, never the child's variants. Regions:
  - the 10 MVA genes: whole gene, 60 kb upstream (the known BUB1B regulatory second alleles
    lie ~44 kb upstream) and 5 kb downstream;
  - the 84 ACMG SF v3.3 genes: MANE Select exons +-20 bp (the screen reads only coding and
    splice variants, so introns would cost ~2 GB of transfer for nothing).
Only frequency fields are kept, so the stored file is small.

Output (git-ignored): refs/gnomad/gnomad_v4.1_joint_slices.vcf.gz (+ .tbi, regions, manifest).
Run: python tools/fetch_gnomad_slices.py [--plan-only]
"""
import argparse
import gzip
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from pipeline.sf_screen import SF_V33  # noqa: E402

MANE_GFF = Path(r"D:\mva-hackathon\refs\ensembl\GRCh38.116.mane_select.gff3.gz")
OUT = REPO / "refs" / "gnomad"
LOG = Path(r"D:\mva-hackathon\logs\gnomad_slices.log")
URL = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{c}.vcf.bgz"
MVA_GENES = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1", "MAD2L1BP", "CEP192", "BUB1"]
ALIASES = {"CENATAC": ["CENATAC", "CCDC84"]}
UPSTREAM, DOWNSTREAM, EXON_PAD = 60_000, 5_000, 20
KEEP = ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AF_grpmax_joint", "grpmax_joint",
        "faf95_joint", "fafmax_faf95_max_joint"]
WSL_TOOLS = "/mnt/d/mva-hackathon/tools/env/bin"
CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_mane():
    genes, tx_gene, exons = {}, {}, defaultdict(list)
    with gzip.open(MANE_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene" and "Name" in a:
                genes[a["Name"]] = (f[0], int(f[3]), int(f[4]), f[6], a.get("ID"))
            elif "Parent" in a and a["Parent"].startswith("gene:") and "ID" in a:
                tx_gene[a["ID"]] = a["Parent"]
            elif f[2] == "exon" and "Parent" in a:
                exons[a["Parent"]].append((int(f[3]), int(f[4])))
    by_gene_id = defaultdict(list)
    for tx, gid in tx_gene.items():
        by_gene_id[gid].extend(exons.get(tx, []))
    return genes, by_gene_id


def build_regions():
    genes, exons_by_gene = read_mane()
    regions, missing, roles = [], [], {}
    for g in MVA_GENES:
        hit = next((n for n in ALIASES.get(g, [g]) if n in genes), None)
        if not hit:
            missing.append(g)
            continue
        chrom, start, end, strand, _ = genes[hit]
        lo, hi = (start - UPSTREAM, end + DOWNSTREAM) if strand == "+" else (start - DOWNSTREAM, end + UPSTREAM)
        regions.append((chrom, max(1, lo), hi))
        roles[hit] = "MVA: gene +60 kb upstream +5 kb downstream"
    for g in sorted(SF_V33):
        if g not in genes:
            missing.append(g)
            continue
        chrom, *_, gid = genes[g]
        ex = exons_by_gene.get(gid, [])
        if not ex:
            missing.append(g)
            continue
        regions += [(chrom, max(1, s - EXON_PAD), e + EXON_PAD) for s, e in ex]
        roles[g] = f"SF v3.3: {len(ex)} MANE exons +-{EXON_PAD} bp"
    merged = defaultdict(list)
    for chrom, s, e in sorted(regions, key=lambda r: (r[0], r[1])):
        if merged[chrom] and s <= merged[chrom][-1][1] + 1:
            merged[chrom][-1][1] = max(merged[chrom][-1][1], e)
        else:
            merged[chrom].append([s, e])
    return merged, roles, missing


def wsl_path(p):
    p = str(p)
    return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plan-only", action="store_true", help="write regions and print totals; fetch nothing")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "idx").mkdir(exist_ok=True)
    (OUT / "parts").mkdir(exist_ok=True)
    merged, roles, missing = build_regions()
    chroms = [c for c in CHROM_ORDER if c in merged]
    total_bp = sum(e - s + 1 for c in chroms for s, e in merged[c])
    for c in chroms:  # 0-based half-open BED with the chr prefix gnomAD uses
        (OUT / "parts" / f"chr{c}.bed").write_text(
            "".join(f"chr{c}\t{s - 1}\t{e}\n" for s, e in merged[c]), encoding="ascii")
    print(f"{len(roles)} genes, {sum(len(merged[c]) for c in chroms)} regions, {total_bp / 1e6:.2f} Mb, "
          f"{len(chroms)} chromosomes; missing from MANE: {missing or 'none'}")
    if args.plan_only:
        return 0

    keep = ",".join(f"INFO/{k}" for k in KEEP)
    parts = []
    lines = ["set -euo pipefail",
             # Every path htslib and bcftools touch is on D:. WSL's own disk image lives on C:.
             "export HOME=/mnt/d/mva-hackathon/tools/home TMPDIR=/mnt/d/DevCache/tmp/wsl",
             'mkdir -p "$TMPDIR"',
             f'cd "{wsl_path(OUT / "idx")}"',
             'LOW="nice -n 19 ionice -c3"']
    for c in chroms:
        part = wsl_path(OUT / "parts" / f"chr{c}.vcf.gz")
        parts.append(part)
        lines.append(
            f'if [ ! -s "{part}" ]; then echo "chr{c}: fetching"; '
            f'$LOW {WSL_TOOLS}/tabix -h -R "{wsl_path(OUT / "parts" / f"chr{c}.bed")}" "{URL.format(c=c)}" '
            f'| $LOW {WSL_TOOLS}/bcftools annotate -x "^{keep}" -Oz -o "{part}.tmp" - '
            f'&& mv "{part}.tmp" "{part}"; else echo "chr{c}: already done"; fi')
    final = wsl_path(OUT / "gnomad_v4.1_joint_slices.vcf.gz")
    lines += [f'{WSL_TOOLS}/bcftools concat -Oz -o "{final}" ' + " ".join(f'"{p}"' for p in parts),
              f'{WSL_TOOLS}/tabix -f -p vcf "{final}"',
              f'echo "records: $({WSL_TOOLS}/bcftools index -n "{final}")"']
    script = OUT / "fetch.sh"
    script.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")

    LOG.parent.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with open(LOG, "a", encoding="utf-8") as log:
        log.write(f"== {started} fetch start\n")
        log.flush()
        proc = subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", wsl_path(script)],
                              stdout=log, stderr=subprocess.STDOUT)
        log.write(f"== exit {proc.returncode}\n")
    final_win = OUT / "gnomad_v4.1_joint_slices.vcf.gz"
    manifest = {
        "approved": "Santanu, live chat 2026-10-04: 'approve gnomAD gene slices', inside the project folder, never C:",
        "source": URL.format(c="{chrom}"), "release": "gnomAD v4.1 joint (exomes + genomes)",
        "fields_kept": KEEP, "genes": roles, "missing_from_mane": missing,
        "regions_bp": total_bp, "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "exit_code": proc.returncode,
        "bytes": final_win.stat().st_size if final_win.exists() else None,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("exit_code", "bytes", "regions_bp", "missing_from_mane")}))
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
