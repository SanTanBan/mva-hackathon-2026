"""Private Kaggle kernel santanubanerjee9/mva-md2pdf: print report HTML to PDF with headless Chromium.

The laptop's Edge printing fails while the laptop is overloaded, so the PDFs are printed here.
Input: the private dataset santanubanerjee9/mva-reports, holding HTML built on the laptop by
tools/render_reports.py's own Markdown-to-HTML step (same CSS as the Edge route). It may hold the
full Track 1 report, which is patient-derived: keep this kernel and that dataset private, and
delete both by 2026-11-23. Output: /kaggle/working/pdf/<name>.pdf, one per HTML file.
"""
import glob
import os
import subprocess
import sys

OUT = "/kaggle/working/pdf"
os.makedirs(OUT, exist_ok=True)
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "playwright"], check=True)
subprocess.run([sys.executable, "-m", "playwright", "install", "--with-deps", "chromium"], check=True,
               stdout=subprocess.DEVNULL)
from playwright.sync_api import sync_playwright  # noqa: E402

files = sorted(h for h in glob.glob("/kaggle/input/**/*.html", recursive=True) if os.path.isfile(h))
print(f"{len(files)} HTML file(s)")
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    for h in files:
        name = os.path.splitext(os.path.basename(h))[0]
        page.goto("file://" + h, wait_until="networkidle")
        page.pdf(path=f"{OUT}/{name}.pdf", format="A4", print_background=True,
                 margin={"top": "18mm", "bottom": "18mm", "left": "16mm", "right": "16mm"})
        print(name, os.path.getsize(f"{OUT}/{name}.pdf"), "bytes")
    browser.close()
