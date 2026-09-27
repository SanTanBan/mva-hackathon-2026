#!/usr/bin/env python3
"""Routine preflight: HF login, gated-dataset access, disk space, and time left.

Prints one JSON object. Reads metadata only and never downloads anything.
"""
import json
import shutil
from datetime import datetime, timezone

REPO = "SageBio/mva-hackathon-2026-data"
CLOSE_UTC = datetime(2026, 10, 24, 23, 59, tzinfo=timezone.utc)
DELETE_BY_UTC = datetime(2026, 11, 23, 23, 59, tzinfo=timezone.utc)


def main():
    now = datetime.now(timezone.utc)
    report = {
        "utc_now": now.isoformat(timespec="minutes"),
        "days_to_close": round((CLOSE_UTC - now).total_seconds() / 86400, 1),
        "days_to_deletion_deadline": round((DELETE_BY_UTC - now).total_seconds() / 86400, 1),
        "disk_free_gb": {d: round(shutil.disk_usage(f"{d}:\\").free / 1e9, 1) for d in ("C", "D")},
    }
    try:
        from huggingface_hub import HfApi
        api = HfApi()
        try:
            report["hf_user"] = api.whoami().get("name")
        except Exception as e:
            report["hf_user"] = None
            report["hf_login_error"] = type(e).__name__
        try:
            api.auth_check(REPO, repo_type="dataset")
            report["dataset_access"] = "granted"
        except Exception as e:
            report["dataset_access"] = "not_granted"
            report["dataset_access_detail"] = f"{type(e).__name__}: {str(e).splitlines()[0][:120]}"
    except ImportError:
        report["dataset_access"] = "unknown"
        report["dataset_access_detail"] = "huggingface_hub is not installed"
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
