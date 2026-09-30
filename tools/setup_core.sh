#!/usr/bin/env bash
# Package R08: install the core toolchain and public references on D: (laptop, VCF-first route).
#
# Approved by Santanu in a live chat on 2026-09-30 (the ~1.5 GB core set). Everything lands under
# D:\mva-hackathon (never on C: or in WSL's home) and is deleted after the competition.
# Run from Windows:  wsl -d Ubuntu -- bash /mnt/d/GitHub_Ideas/For_Profit/Child-with-rare-Disease/tools/setup_core.sh
# Steps are idempotent: a finished step is skipped on rerun. Heavy steps run under nice/ionice.
set -euo pipefail

BASE=/mnt/d/mva-hackathon
TOOLS=$BASE/tools
REFS=$BASE/refs
LOG=$BASE/logs/r08_setup_core.log
export HOME=$TOOLS/home XDG_CACHE_HOME=$TOOLS/cache MAMBA_ROOT_PREFIX=$TOOLS/mamba CONDA_PKGS_DIRS=$TOOLS/mamba/pkgs
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$MAMBA_ROOT_PREFIX" "$TOOLS/bin" "$REFS/grch38" "$REFS/clinvar" "$REFS/ensembl" "$(dirname "$LOG")"
exec >> "$LOG" 2>&1
LOW="nice -n 19 ionice -c3"
echo "== $(date -u +%FT%TZ) setup_core start"

# 1. micromamba: a single static binary from the official release channel.
if [ ! -x "$TOOLS/bin/micromamba" ]; then
  curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$TOOLS" bin/micromamba
fi
"$TOOLS/bin/micromamba" --version | sed 's/^/micromamba /'

# 2. Genomics toolchain from conda-forge and bioconda. --always-copy avoids hard links on the D: mount.
if [ ! -x "$TOOLS/env/bin/bcftools" ]; then
  $LOW "$TOOLS/bin/micromamba" create -y -q -p "$TOOLS/env" --always-copy -c conda-forge -c bioconda \
    bcftools samtools htslib minimap2
fi
for t in bcftools samtools bgzip tabix minimap2; do
  printf '%s: ' "$t"; "$TOOLS/env/bin/$t" --version 2>&1 | head -1
done

# 3. GRCh38 no-alt analysis set with the hs38d1 decoys (NCBI), contigs renamed without "chr" to
#    match the dataset VCF, then indexed.
FA=$REFS/grch38/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa
if [ ! -s "$FA.fai" ]; then
  URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
  curl -fsSL --retry 3 -o "$REFS/grch38/ref.fna.gz" "$URL"
  $LOW bash -c "gzip -dc '$REFS/grch38/ref.fna.gz' | sed 's/^>chr/>/' > '$FA.tmp'"
  mv "$FA.tmp" "$FA"
  rm -f "$REFS/grch38/ref.fna.gz"
  $LOW "$TOOLS/env/bin/samtools" faidx "$FA"
fi
echo "grch38: $(wc -l < "$FA.fai") contigs, $(du -h "$FA" | cut -f1)"

# 4. ClinVar GRCh38 VCF and its index (NCBI). The release date is kept from the file header.
if [ ! -s "$REFS/clinvar/clinvar.vcf.gz.tbi" ]; then
  curl -fsSL --retry 3 -o "$REFS/clinvar/clinvar.vcf.gz" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  curl -fsSL --retry 3 -o "$REFS/clinvar/clinvar.vcf.gz.tbi" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi
echo "clinvar: $("$TOOLS/env/bin/bcftools" view -h "$REFS/clinvar/clinvar.vcf.gz" | grep -m1 '^##fileDate')"

# 5. Ensembl GRCh38 GFF3, current release.
if ! ls "$REFS"/ensembl/Homo_sapiens.GRCh38.*.gff3.gz >/dev/null 2>&1; then
  GFF=$(curl -fsSL https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/ | grep -o 'Homo_sapiens\.GRCh38\.[0-9]*\.gff3\.gz' | head -1)
  curl -fsSL --retry 3 -o "$REFS/ensembl/$GFF" "https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/$GFF"
fi
echo "ensembl: $(ls "$REFS"/ensembl/)"

du -sh "$TOOLS" "$REFS" | sed 's/^/size /'
echo "== $(date -u +%FT%TZ) setup_core done"
