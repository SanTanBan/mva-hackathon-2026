#!/usr/bin/env bash
# Package D02 annotation step, run in WSL: normalise, MANE Select consequences, ClinVar fields.
#   1. bcftools norm: split multi-allelics and left-align indels against the no-chr GRCh38
#      analysis set (the VCF has no 'chr' prefix; the reference, GFF3 and ClinVar match it).
#   2. bcftools csq: consequences on MANE Select transcripts only (pipeline/mane_gff.py),
#      written to INFO/BCSQ.
#   3. bcftools annotate: CLNSIG, CLNREVSTAT, ALLELEID and the ClinVar population frequencies
#      (AF_ESP, AF_EXAC, AF_TGP), matched on position and alleles.
#   4. When refs/gnomad holds the gnomAD v4.1 slices and their covered-intervals BED
#      (tools/fetch_gnomad_slices.py, package R10): AF_joint, nhomalt_joint and
#      fafmax_faf95_max_joint, matched on chrom, pos, ref and alt, plus the flag
#      INFO/GNOMAD_COVERED for records inside a fetched interval.
# pipeline/annotate_table.py then turns the output into the table pipeline/comphet.py reads.
#
# Usage (from Windows):
#   wsl -d Ubuntu -- bash /mnt/d/GitHub_Ideas/For_Profit/Child-with-rare-Disease/tools/annotate_vcf.sh \
#       /mnt/d/mva-hackathon/data/IN.vcf.gz /mnt/d/mva-hackathon/work/annotated.vcf.gz
# Output must be outside the repository (patient-derived). Everything runs at low priority;
# temporary files go next to the output on D:, never to WSL's home or /tmp (both on C:).
set -euo pipefail
IN=${1:?input VCF}; OUT=${2:?output .vcf.gz}
case "$(realpath -m "$OUT")" in
  /mnt/d/GitHub_Ideas/*) echo "refusing to write inside the repository: $OUT" >&2; exit 2;;
esac
TOOLS=/mnt/d/mva-hackathon/tools/env/bin
REFS=/mnt/d/mva-hackathon/refs
FA=${FA:-$REFS/grch38/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa}
GFF=${GFF:-$REFS/ensembl/GRCh38.116.mane_select.gff3.gz}
CLINVAR=${CLINVAR:-$REFS/clinvar/clinvar.vcf.gz}
GNOMAD_DIR=$(dirname "$(realpath "$0")")/../refs/gnomad
GNOMAD=${GNOMAD:-$GNOMAD_DIR/gnomad_v4.1_joint_slices.vcf.gz}
[ -f "$GNOMAD" ] || GNOMAD=$GNOMAD_DIR/gnomad_v4.1_joint_slices.vcf  # compression deferred
GNOMAD_BED=${GNOMAD_BED:-$GNOMAD_DIR/gnomad_v4.1_covered.bed}
LOW="nice -n 19 ionice -c3"
export TMPDIR; TMPDIR=$(dirname "$(realpath -m "$OUT")")/tmp.annotate.$$
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT

$LOW "$TOOLS/bcftools" norm -m-any -f "$FA" --check-ref w -Ou "$IN" \
  | $LOW "$TOOLS/bcftools" csq -f "$FA" -g "$GFF" -p a --ncsq 64 -Ob -o "$TMPDIR/csq.bcf" 2> "$TMPDIR/csq.log"
grep -E "^(Indexed|Calling)" "$TMPDIR/csq.log" || true
# annotate -a reads both files by region, so its input must be an indexed file, not a pipe.
"$TOOLS/bcftools" index "$TMPDIR/csq.bcf"
if [ -f "$GNOMAD" ] && [ -f "$GNOMAD_BED" ]; then
  $LOW "$TOOLS/bcftools" annotate -a "$CLINVAR" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/ALLELEID,INFO/AF_ESP,INFO/AF_EXAC,INFO/AF_TGP \
    -Ob -o "$TMPDIR/clinvar.bcf" "$TMPDIR/csq.bcf"
  "$TOOLS/bcftools" index "$TMPDIR/clinvar.bcf"
  # The slices use chr contigs; the VCF does not. Rename a copy, then match on chrom, pos, ref
  # and alt (the slices are already split and left-aligned by gnomAD).
  for c in $(seq 1 22) X Y; do printf 'chr%s\t%s\n' "$c" "$c"; done > "$TMPDIR/chrmap.txt"
  $LOW "$TOOLS/bcftools" annotate --rename-chrs "$TMPDIR/chrmap.txt" -Ob -o "$TMPDIR/gnomad.bcf" "$GNOMAD"
  "$TOOLS/bcftools" index "$TMPDIR/gnomad.bcf"
  "$TOOLS/bgzip" -c "$GNOMAD_BED" > "$TMPDIR/covered.bed.gz"
  "$TOOLS/tabix" -p bed "$TMPDIR/covered.bed.gz"
  $LOW "$TOOLS/bcftools" annotate -a "$TMPDIR/gnomad.bcf" \
    -c INFO/AF_joint,INFO/nhomalt_joint,INFO/fafmax_faf95_max_joint -Ou "$TMPDIR/clinvar.bcf" \
    | $LOW "$TOOLS/bcftools" annotate -a "$TMPDIR/covered.bed.gz" -c CHROM,FROM,TO -m +GNOMAD_COVERED \
        -Oz -o "$OUT"
else
  echo "gnomAD slices not found ($GNOMAD, $GNOMAD_BED): ClinVar frequencies only" >&2
  $LOW "$TOOLS/bcftools" annotate -a "$CLINVAR" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/ALLELEID,INFO/AF_ESP,INFO/AF_EXAC,INFO/AF_TGP \
    -Oz -o "$OUT" "$TMPDIR/csq.bcf"
fi
"$TOOLS/bcftools" index -t "$OUT"
echo "records: $("$TOOLS/bcftools" view -H "$OUT" | wc -l)"
echo "with BCSQ: $("$TOOLS/bcftools" view -H -i 'INFO/BCSQ!="."' "$OUT" | wc -l)"
echo "in ClinVar: $("$TOOLS/bcftools" view -H -i 'INFO/ALLELEID!="."' "$OUT" | wc -l)"
if "$TOOLS/bcftools" view -h "$OUT" | grep -q 'ID=GNOMAD_COVERED'; then
  echo "in gnomAD slices: $("$TOOLS/bcftools" view -H -i 'INFO/GNOMAD_COVERED=1' "$OUT" | wc -l)"
  echo "listed by gnomAD: $("$TOOLS/bcftools" view -H -i 'INFO/AF_joint!="."' "$OUT" | wc -l)"
fi
"$TOOLS/bcftools" --version | head -1
