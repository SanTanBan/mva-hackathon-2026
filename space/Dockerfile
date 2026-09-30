FROM python:3.11-slim

# Genomics toolchain from Debian: bcftools, samtools, htslib (tabix, bgzip) and minimap2; poppler for PDF text.
RUN apt-get update \
    && apt-get install -y --no-install-recommends bcftools samtools tabix minimap2 poppler-utils curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH PYTHONUNBUFFERED=1
WORKDIR /home/user/app

COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt
COPY --chown=user . .

EXPOSE 7860
CMD ["python", "app.py"]
