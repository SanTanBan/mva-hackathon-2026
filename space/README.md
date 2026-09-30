---
title: MVA private runner
emoji: 🧬
colorFrom: gray
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
---

# MVA Hackathon private runner (package R08)

This is a private Space. It reads the gated MVA Hackathon dataset in place, using an `HF_TOKEN` secret that the owner adds himself, so the child's genome never leaves Hugging Face. It returns only small result tables to the owner. Everything is deleted by 2026-11-23, as the Hackathon rules require.

The source lives in the project repository under `space/`. The pipeline code is uploaded alongside it from `pipeline/`.
