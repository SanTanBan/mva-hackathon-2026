#!/usr/bin/env python3
"""Run the Track 1 fast path in the private Kaggle kernel santanubanerjee9/mva-d02 and collect it.

  1. push kaggle/mva_d02 (unless --collect-only);
  2. poll its status every 90 s;
  3. download only the small outputs (d02/*.csv, *.tsv, summary.json and the log) to
     D:\\mva-hackathon\\work\\d02\\<date>, never inside the repository;
  4. validate the CSV with tools/validate_track1_csv.py.
Patient-derived output: nothing from it is printed except counts and validator messages.

Run: python tools/kaggle_d02.py [--collect-only] [--max-minutes 150]
"""
import argparse
import os
import re
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.common import assert_outside_repo  # noqa: E402
from tools.validate_track1_csv import validate  # noqa: E402

KERNEL = "santanubanerjee9/mva-d02"
KTMP = Path(r"D:\DevCache\tmp\claude\Child-with-rare-Disease\ktmp")
HF_USER = "SantanuBanerjee"


def kaggle(*args):
    KTMP.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TEMP=str(KTMP), TMP=str(KTMP))
    r = subprocess.run(["kaggle", *args], env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")  # progress bars are not cp1252
    out = "\n".join(ln for ln in (r.stdout + r.stderr).splitlines()
                    if "warn" not in ln.lower() and "%|" not in ln)
    return r.returncode, out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--collect-only", action="store_true")
    ap.add_argument("--max-minutes", type=int, default=150)
    args = ap.parse_args()
    if not args.collect_only:
        rc, out = kaggle("kernels", "push", "-p", str(ROOT / "kaggle" / "mva_d02"))
        print(out.strip())
        if rc:
            return 1
        time.sleep(30)
    deadline = time.time() + args.max_minutes * 60
    while True:
        _, out = kaggle("kernels", "status", KERNEL)
        status = (re.findall(r"KernelWorkerStatus\.([A-Z_]+)", out) or ["UNKNOWN"])[0]
        print(time.strftime("%H:%M:%S"), status, flush=True)
        if status not in ("RUNNING", "QUEUED", "UNKNOWN"):
            break
        if time.time() > deadline:
            print("still running; collect later with --collect-only")
            return 2
        time.sleep(90)
    dest = assert_outside_repo(Path(r"D:\mva-hackathon\work\d02") / date.today().isoformat())
    dest.mkdir(parents=True, exist_ok=True)
    rc, out = kaggle("kernels", "output", KERNEL, "-p", str(dest), "-o", "-q",
                     "--file-pattern", r".*(\.csv|\.tsv|summary\.json|\.log)$")
    print(out.strip()[-500:])
    files = sorted(p.relative_to(dest).as_posix() for p in dest.rglob("*") if p.is_file())
    print("collected:", ", ".join(files) or "nothing")
    csvs = list(dest.rglob(f"{HF_USER}_track1.csv"))
    if status != "COMPLETE" or not csvs:
        print(f"kernel status {status}; read the log in {dest}")
        return 1
    errors, warnings = validate(csvs[0], hf_user=HF_USER)
    print("validator errors:", errors or "none")
    print("validator warnings:", warnings or "none")
    print("csv:", csvs[0])
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
