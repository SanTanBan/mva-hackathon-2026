#!/usr/bin/env python3
"""Package D02 fast path, in one command, once D00 and D01 are done.

  1. tools/annotate_vcf.sh in WSL: normalise, MANE Select consequences, ClinVar fields
     (about 5 million records; expect well over 20 minutes on this laptop, so run it in the
     background with output to D:\\mva-hackathon\\logs).
  2. pipeline/annotate_table.py: the per-gene table.
  3. pipeline/track1_rank.py: the ranked Track 1 CSV.
  4. tools/validate_track1_csv.py on the CSV.
Everything is written to D:\\mva-hackathon\\work\\d02. Staging into
D:\\mva-hackathon\\submissions\\track1\\<date> (with the report and SUBMIT.md) is a separate,
reviewed step, and submitting needs Santanu's live yes.

Run: python tools/run_d02.py --proband-id ID [--vcf PATH] [--panel-genes FILE] [--skip-annotate]
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import comphet, track1_rank  # noqa: E402
from pipeline.common import assert_outside_repo  # noqa: E402
from tools.validate_track1_csv import validate  # noqa: E402

WORK = Path(r"D:\mva-hackathon\work")
OUT = WORK / "d02"
GFF = Path(r"D:\mva-hackathon\refs\ensembl\GRCh38.116.mane_select.gff3.gz")
PANEL_GENES = Path(r"D:\mva-hackathon\refs\panelapp\panel_genes_mva.txt")  # tools/build_panel_genes.py
HF_USER = "SantanuBanerjee"


def wsl_path(p):
    p = Path(p).resolve()
    return "/mnt/" + p.drive[0].lower() + p.as_posix()[2:]


def find_vcf():
    inventory = json.loads((WORK / "d00_inventory.json").read_text(encoding="utf-8"))
    vcfs = [f["local"] for f in inventory["files"] if f.get("local") and f["kind"] == "vcf"]
    if len(vcfs) != 1:
        raise SystemExit(f"expected one VCF in the D00 inventory, found {len(vcfs)}; pass --vcf")
    return Path(vcfs[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proband-id", required=True)
    ap.add_argument("--vcf")
    ap.add_argument("--panel-genes", default=str(PANEL_GENES) if PANEL_GENES.exists() else None,
                    help="one gene symbol per line (default: the pre-registered PanelApp file, if built)")
    ap.add_argument("--skip-annotate", action="store_true", help="reuse an existing annotated VCF")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    assert_outside_repo(OUT)
    annotated = OUT / "annotated.vcf.gz"
    if not args.skip_annotate:
        vcf = Path(args.vcf) if args.vcf else find_vcf()
        script = wsl_path(Path(__file__).with_name("annotate_vcf.sh"))
        r = subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", script, wsl_path(vcf), wsl_path(annotated)],
                           env=dict(os.environ, MSYS_NO_PATHCONV="1"))
        if r.returncode:
            raise SystemExit(f"annotate_vcf.sh failed ({r.returncode})")
    table = OUT / "annotated_table.tsv"
    subprocess.run([sys.executable, "-m", "pipeline.annotate_table", str(annotated), str(GFF), str(table)],
                   check=True, cwd=Path(__file__).resolve().parents[1])
    panel = []
    if args.panel_genes:
        panel = [ln.strip() for ln in open(args.panel_genes, encoding="utf-8") if ln.strip() and not ln.startswith("#")]
        print(f"panel genes: {len(panel)} from {args.panel_genes}")
    else:
        print("panel genes: none (run tools/build_panel_genes.py first)")
    ranked = track1_rank.score_candidates(comphet.load(table), panel)
    csv_path = OUT / f"{HF_USER}_track1.csv"
    track1_rank.write_csv(csv_path, args.proband_id, ranked)
    errors, warnings = validate(csv_path, hf_user=HF_USER)
    print(f"{len(ranked)} rows -> {csv_path}")
    print("validator errors:", errors or "none")
    print("validator warnings:", warnings or "none")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
