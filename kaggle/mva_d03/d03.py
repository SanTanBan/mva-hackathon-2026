"""Private Kaggle kernel santanubanerjee9/mva-d03: low-pass depth for mosaic aneuploidy (package D03).

PATIENT DATA (FASTQ). Runs only as a private mva-* kernel on santanubanerjee9 (CLAUDE.md, approvals of
2026-10-14 and 2026-10-18: FASTQ steps run only inside private Kaggle/HF compute, never on the laptop).

Stage 1 (this version): check that the Kaggle secret HF_TOKEN is attached and can read the gated
dataset, and list the FASTQ file sizes. No read data is downloaded.

The token is Santanu's own fine-grained, read-only Hugging Face token, added by him as a Kaggle secret.
It is never printed, written to output, or passed anywhere except the Hugging Face API.
"""
import json
import urllib.request

REPO = "SageBio/mva-hackathon-2026-data"
out = {"stage": "access check"}
try:
    from kaggle_secrets import UserSecretsClient
    token = UserSecretsClient().get_secret("HF_TOKEN")
    out["secret_attached"] = bool(token)
except Exception as e:  # secret missing or not attached to this notebook
    token, out["secret_attached"], out["secret_error"] = None, False, type(e).__name__

if token:
    req = urllib.request.Request(f"https://huggingface.co/api/datasets/{REPO}/tree/main?recursive=true",
                                 headers={"Authorization": f"Bearer {token}"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            tree = json.loads(r.read().decode())
        files = [x for x in tree if x.get("type") == "file"]
        out["access"] = "ok"
        out["fastq_gb"] = {x["path"].rsplit("/", 1)[-1][-24:]: round(x.get("size", 0) / 1e9, 2)
                           for x in files if ".fastq" in x["path"]}
        out["fastq_total_gb"] = round(sum(x.get("size", 0) for x in files if ".fastq" in x["path"]) / 1e9, 1)
    except Exception as e:
        out["access"] = f"failed: {type(e).__name__} {getattr(e, 'code', '')}"
print(json.dumps(out, indent=1))
json.dump(out, open("/kaggle/working/d03_access.json", "w"), indent=1)
