"""Private runner for the MVA Hackathon gated data (package R08).

Runs inside a private Hugging Face Space owned by SantanuBanerjee. It reads the gated dataset
SageBio/mva-hackathon-2026-data in place, using the HF_TOKEN secret that Santanu adds himself, so the
child's genome never leaves Hugging Face (CLAUDE.md, amended 2026-09-30). Nothing here sends patient
data to any other service. Work files live on the Space's ephemeral disk under /tmp and disappear when
the Space restarts; everything is deleted by 2026-11-23 in any case.

Actions (each is also exposed as a Gradio API endpoint for the routine to call):
- self_test: tool versions, the pipeline's synthetic tests, and a synthetic bcftools round trip.
- check_access: whether HF_TOKEN can read the gated dataset; file names and sizes only.
- prepare_reference: GRCh38 no-alt analysis set with hs38d1 decoys from NCBI, contigs renamed without
  the chr prefix to match the dataset VCF, then indexed with samtools faidx.
- d01_intake: inventory with checksums, VCF header summary, and HPO terms from the phenotype
  document. Refuses to run unless access is granted.
"""
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
import unittest
from pathlib import Path

import gradio as gr
from huggingface_hub import HfApi, hf_hub_download

DATASET = "SageBio/mva-hackathon-2026-data"
WORK, REFS = Path("/tmp/work"), Path("/tmp/refs")
REF_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
           "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz")
REF_FA = REFS / "GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa"
TOOLS = [["bcftools", "--version"], ["samtools", "--version"], ["bgzip", "--version"], ["minimap2", "--version"]]
SELF_TESTS = ["pipeline.tests.test_pipeline", "pipeline.tests.test_allele_mechanism",
              "pipeline.tests.test_sf_screen", "pipeline.tests.test_reads"]
HPO = re.compile(r"HP:\d{7}")


def sh(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    return r.returncode, (r.stdout + r.stderr).strip()


def token():
    return os.environ.get("HF_TOKEN", "").strip()


def access_granted():
    if not token():
        return False, "The HF_TOKEN secret is not set. Santanu adds it in the Space settings."
    try:
        HfApi(token=token()).auth_check(DATASET, repo_type="dataset")
        return True, "granted"
    except Exception as e:  # GatedRepoError, RepositoryNotFoundError, HTTP errors
        return False, f"not granted: {type(e).__name__}"


def self_test():
    lines = [f"python {sys.version.split()[0]}"]
    for cmd in TOOLS:
        rc, out = sh(cmd)
        lines.append(f"{cmd[0]}: {out.splitlines()[0] if rc == 0 and out else 'MISSING'}")
    buf = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromNames(SELF_TESTS)
    res = unittest.TextTestRunner(stream=buf, verbosity=1).run(suite)
    lines.append(f"pipeline tests: {res.testsRun} run, {len(res.failures)} failures, {len(res.errors)} errors, "
                 f"{len(res.skipped)} skipped")
    if res.failures or res.errors:
        lines.append(buf.getvalue()[-3000:])
    # Synthetic round trip: a two-record VCF on a no-chr contig, compressed, indexed and queried.
    d = WORK / "selftest"
    d.mkdir(parents=True, exist_ok=True)
    vcf = d / "synthetic.vcf"
    rec = lambda pos, ref, alt, gt: "\t".join(["15", str(pos), ".", ref, alt, "50", "PASS", ".", "GT", gt])
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=15,length=1000000>\n"
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="synthetic">\n'
                   + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SYNTH"]) + "\n"
                   + rec(1000, "A", "G", "0/1") + "\n" + rec(2000, "C", "T", "1/1") + "\n")
    ok = sh(["bgzip", "-f", str(vcf)])[0] == 0 and sh(["tabix", "-f", "-p", "vcf", f"{vcf}.gz"])[0] == 0
    rc, out = sh(["bcftools", "view", "-H", "-r", "15:1500-2500", f"{vcf}.gz"])
    lines.append("bcftools round trip: " + ("ok" if ok and rc == 0 and out.count("\n") == 0 and "\t2000\t" in out else f"FAILED {out[:200]}"))
    shutil.rmtree(d, ignore_errors=True)
    granted, why = access_granted()
    lines.append(f"dataset access: {why}")
    lines.append(f"reference ready: {REF_FA.exists() and Path(str(REF_FA) + '.fai').exists()}")
    rc, out = sh(["df", "-h", "/tmp"])
    lines.append(out.splitlines()[-1] if rc == 0 else "df unavailable")
    return "\n".join(lines)


def check_access():
    granted, why = access_granted()
    if not granted:
        return why
    api = HfApi(token=token())
    rows = []
    for f in api.list_repo_tree(DATASET, repo_type="dataset", recursive=True):
        size = getattr(f, "size", None)
        if size is not None:
            rows.append(f"{f.path}\t{size / 1e9:.2f} GB")
    return f"access granted; {len(rows)} files\n" + "\n".join(rows)


def prepare_reference():
    if REF_FA.exists() and Path(str(REF_FA) + ".fai").exists():
        return f"already prepared: {REF_FA}"
    REFS.mkdir(parents=True, exist_ok=True)
    gz = REFS / "ref.fna.gz"
    t0 = time.time()
    rc, out = sh(["curl", "-fsSL", "--retry", "3", "-o", str(gz), REF_URL])
    if rc != 0:
        return f"download failed: {out[:300]}"
    with gzip.open(gz, "rt") as src, open(REF_FA, "w") as dst:
        for line in src:
            dst.write(">" + line[4:] if line.startswith(">chr") else line)
    gz.unlink()
    rc, out = sh(["samtools", "faidx", str(REF_FA)])
    if rc != 0:
        return f"faidx failed: {out[:300]}"
    n = sum(1 for _ in open(str(REF_FA) + ".fai"))
    return f"reference ready in {time.time() - t0:.0f} s: {n} contigs, {REF_FA.stat().st_size / 1e9:.2f} GB"


def d01_intake():
    """Inventory, VCF header summary and HPO terms. Returns JSON for D:\\mva-hackathon\\work only."""
    granted, why = access_granted()
    if not granted:
        return json.dumps({"status": "refused", "reason": why})
    api = HfApi(token=token())
    files = []
    for f in api.list_repo_tree(DATASET, repo_type="dataset", recursive=True):
        if getattr(f, "size", None) is None:
            continue
        lfs = getattr(f, "lfs", None)
        files.append({"path": f.path, "bytes": f.size, "sha256": getattr(lfs, "sha256", None) if lfs else None,
                      "blob_id": getattr(f, "blob_id", None)})
    out = {"status": "ok", "files": files, "vcf": [], "hpo": [], "phenotype_files": []}
    WORK.mkdir(parents=True, exist_ok=True)
    for f in files:
        p = f["path"]
        if p.endswith((".vcf.gz", ".vcf")):
            local = hf_hub_download(DATASET, p, repo_type="dataset", token=token(), cache_dir=str(WORK / "hf"))
            rc, hdr = sh(["bcftools", "view", "-h", local])
            contigs = re.findall(r"##contig=<ID=([^,>]+)", hdr)
            out["vcf"].append({
                "path": p,
                "samples": hdr.splitlines()[-1].split("\t")[9:] if rc == 0 else [],
                "reference": re.findall(r"##reference=(.*)", hdr),
                "n_contigs": len(contigs), "first_contigs": contigs[:5],
                "chr_prefix": any(c.startswith("chr") for c in contigs),
                "info": re.findall(r"##INFO=<ID=([^,]+)", hdr), "format": re.findall(r"##FORMAT=<ID=([^,]+)", hdr),
                "sources": re.findall(r"##source=(.*)", hdr)[:5],
                "n_records": int(sh(["bash", "-c", f"bcftools view -H '{local}' | wc -l"])[1] or 0),
            })
        elif not p.endswith((".fastq.gz", ".fq.gz", ".tbi", ".csi", ".bam", ".cram", ".bai", ".crai", ".gitattributes")):
            local = hf_hub_download(DATASET, p, repo_type="dataset", token=token(), cache_dir=str(WORK / "hf"))
            text = ""
            if p.lower().endswith(".pdf"):
                text = sh(["pdftotext", "-layout", local, "-"])[1]
            else:
                try:
                    text = Path(local).read_text(encoding="utf-8", errors="replace")
                except Exception:
                    text = ""
            out["phenotype_files"].append({"path": p, "chars": len(text)})
            out["hpo"] += sorted(set(HPO.findall(text)) - set(out["hpo"]))
    return json.dumps(out, indent=1)


with gr.Blocks(title="MVA private runner") as demo:
    gr.Markdown("### MVA Hackathon private runner\nPrivate Space. Reads the gated dataset in place with the HF_TOKEN "
                "secret. No patient data leaves Hugging Face except small results fetched by the owner.")
    out = gr.Textbox(label="Output", lines=24)
    with gr.Row():
        gr.Button("Self-test").click(self_test, outputs=out, api_name="self_test")
        gr.Button("Check access").click(check_access, outputs=out, api_name="check_access")
        gr.Button("Prepare reference").click(prepare_reference, outputs=out, api_name="prepare_reference")
        gr.Button("D01 intake").click(d01_intake, outputs=out, api_name="d01_intake")

if __name__ == "__main__":
    demo.queue(default_concurrency_limit=1).launch(server_name="0.0.0.0", server_port=7860)
