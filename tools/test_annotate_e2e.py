#!/usr/bin/env python3
"""End-to-end check of the D02 annotation path on a synthetic VCF, with the real toolchain.

The VCF is built at run time from public ClinVar BUB1B variants (cited by position in ClinVar,
not from the child) with invented genotypes, plus two invented SNVs (one upstream, one deep intronic). It runs
tools/annotate_vcf.sh in WSL, then pipeline/annotate_table.py and pipeline/comphet.py, and checks:
  - the multi-allelic record is split and both alleles are annotated;
  - MANE Select consequences: stop_gained, frameshift, splice_region, deep intron, synonymous,
    upstream;
  - ClinVar CLNSIG is attached to the known pathogenic alleles;
  - comphet pairs the two read-backed heterozygous alleles as 'trans';
  - pipeline/track1_rank.py ranks that pair first and writes a CSV that passes the validator.
Everything is written under D:\\DevCache\\tmp\\claude\\mva and deleted afterwards.
It is skipped (exit 0) when WSL or the toolchain is missing.

Run: python tools/test_annotate_e2e.py
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import annotate_table, comphet, track1_rank  # noqa: E402
from tools.validate_track1_csv import validate  # noqa: E402
from pipeline.common import open_text  # noqa: E402

TMP = Path(r"D:\DevCache\tmp\claude\mva\e2e")
GFF = Path(r"D:\mva-hackathon\refs\ensembl\GRCh38.116.mane_select.gff3.gz")
TOOLS = Path(r"D:\mva-hackathon\tools\env\bin\bcftools")
REPO_WSL = "/mnt/d/GitHub_Ideas/For_Profit/Child-with-rare-Disease"
TMP_WSL = "/mnt/d/DevCache/tmp/claude/mva/e2e"


def t(*f):
    return "\t".join(str(x) for x in f)


def wsl(cmd):
    env = dict(os.environ, MSYS_NO_PATHCONV="1")
    return subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", "-c", cmd], capture_output=True,
                          text=True, env=env, timeout=1800)


def ref_base(pos):
    fa = "/mnt/d/mva-hackathon/refs/grch38/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa"
    r = wsl(f"/mnt/d/mva-hackathon/tools/env/bin/samtools faidx {fa} 15:{pos}-{pos} | tail -1")
    return r.stdout.strip().upper()


def build_vcf(path):
    up = 40158000  # 2,984 bp upstream of the BUB1B MANE Select gene start (40,160,984, + strand)
    deep = 40168000  # deep in a BUB1B intron
    b, d = ref_base(up), ref_base(deep)
    alt, dalt = ("A" if b != "A" else "G"), ("A" if d != "A" else "G")
    fmt = "GT:AD:DP:GQ:PGT:PID"
    pid = "40165186_C_T"
    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38_no_alt_analysis_set_no_chr (synthetic test)",
        "##contig=<ID=15,length=101991189>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype">',
        '##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID">',
        t("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SYNTH"),
    ]
    recs = [
        t(15, up, ".", b, alt, 60, "PASS", ".", fmt, "0/1:15,14:29:99:.:."),
        t(15, 40161226, ".", "G", "A", 60, "PASS", ".", fmt, "0/1:15,15:30:99:.:."),
        t(15, 40165117, ".", "C", "T,A", 60, "PASS", ".", fmt, "1/2:0,12,13:25:99:.:."),
        t(15, 40165186, ".", "C", "T", 60, "PASS", ".", fmt, f"0/1:16,14:30:99:0|1:{pid}"),
        t(15, 40165201, ".", "G", "A", 60, "PASS", ".", fmt, "0/1:13,12:25:99:.:."),
        t(15, deep, ".", d, dalt, 60, "PASS", ".", fmt, "0/1:12,12:24:99:.:."),
        t(15, 40170578, ".", "AGG", "A", 60, "PASS", ".", fmt, f"0/1:14,13:27:99:1|0:{pid}"),
    ]
    path.write_text("\n".join(header + recs) + "\n", encoding="utf-8", newline="\n")


def main():
    if not (TOOLS.exists() and GFF.exists() and shutil.which("wsl")):
        print("SKIP: toolchain, MANE GFF3 or WSL missing")
        return 0
    TMP.mkdir(parents=True, exist_ok=True)
    try:
        build_vcf(TMP / "synth.vcf")
        r = wsl(f"cd {TMP_WSL} && /mnt/d/mva-hackathon/tools/env/bin/bgzip -f synth.vcf && "
                f"bash {REPO_WSL}/tools/annotate_vcf.sh {TMP_WSL}/synth.vcf.gz {TMP_WSL}/annotated.vcf.gz")
        print(r.stdout.strip())
        if r.returncode:
            print(r.stderr)
            return 1
        with open_text(GFF) as fh:
            genes = annotate_table.GeneIndex(fh)
        with open_text(TMP / "annotated.vcf.gz") as fh:
            table = [row for row in annotate_table.rows(fh, genes) if row["gene"] == "BUB1B"]
        by_pos = {(row["pos"], row["alt"]): row for row in table}
        for row in table:
            print(t(row["pos"], row["ref"][:6], row["alt"][:6], row["consequence"], row["impact"],
                    row["aa_change"], row["clnsig"], row["pop_af"]))
        checks = {
            "split multi-allelic": ("40165117", "T") in by_pos and ("40165117", "A") in by_pos,
            "stop_gained": by_pos[("40165186", "T")]["consequence"] == "stop_gained",
            "frameshift": by_pos[("40170578", "A")]["consequence"].startswith("frameshift"),
            "splice_region": by_pos[("40165201", "A")]["consequence"] == "splice_region",
            "deep intron": [r["consequence"] for r in table if r["pos"] == "40168000"] in (["intron"], ["intron_variant"]),
            "synonymous": by_pos[("40161226", "A")]["consequence"] == "synonymous",
            "upstream": any(r["consequence"] == "upstream_gene_variant" for r in table),
            "clinvar attached": "athogenic" in by_pos[("40165186", "T")]["clnsig"],
            "chr prefix": all(r["chrom"] == "chr15" for r in table),
        }
        cands = comphet.candidates([{**r, "gene": r["gene"]} for r in table])
        pair = [c for c in cands if c["variant_1"] == "chr15:40165186:C>T" and c["variant_2"] == "chr15:40170578:AGG>A"]
        checks["comphet trans"] = bool(pair) and pair[0]["phase"] == "trans"
        ranked = track1_rank.score_candidates(table)
        csv_path = TMP / "SantanuBanerjee_track1_synthetic.csv"
        track1_rank.write_csv(csv_path, "SYNTHETIC-1", ranked)
        print(csv_path.read_text(encoding="utf-8"))
        errors, warnings = validate(csv_path, hf_user="SantanuBanerjee")
        print("validator:", errors or "no errors", warnings or "")
        checks["csv valid"] = not errors
        checks["trans pathogenic pair ranks first"] = (ranked[0][1]["pos"], ranked[0][2]["pos"]) == ("40165186", "40170578")
        for name, ok in checks.items():
            print(("PASS " if ok else "FAIL ") + name)
        return 0 if all(checks.values()) else 1
    finally:
        shutil.rmtree(TMP, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
