#!/usr/bin/env bash
# Package D02 annotation step, run in WSL: normalise, MANE Select consequences, ClinVar fields.
#   1. bcftools norm: split multi-allelics and left-align indels against the no-chr GRCh38
#      analysis set (the VCF has no 'chr' prefix; the reference, GFF3 and ClinVar match it).
#   2. bcftools csq: consequences on MANE Select transcripts only (pipeline/mane_gff.py),
#      written to INFO/BCSQ.
#   3. bcftools annotate: CLNSIG, CLNREVSTAT, ALLELEID and the ClinVar population frequencies
#      (AF_ESP, AF_EXAC, AF_TGP), matched on position and alleles.
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
LOW="nice -n 19 ionice -c3"
export TMPDIR; TMPDIR=$(dirname "$(realpath -m "$OUT")")/tmp.annotate.$$
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT

$LOW "$TOOLS/bcftools" norm -m-any -f "$FA" --check-ref w -Ou "$IN" \
  | $LOW "$TOOLS/bcftools" csq -f "$FA" -g "$GFF" -p a --ncsq 64 -Ob -o "$TMPDIR/csq.bcf" 2> "$TMPDIR/csq.log"
grep -E "^(Indexed|Calling)" "$TMPDIR/csq.log" || true
# annotate -a reads both files by region, so its input must be an indexed file, not a pipe.
"$TOOLS/bcftools" index "$TMPDIR/csq.bcf"
$LOW "$TOOLS/bcftools" annotate -a "$CLINVAR" \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/ALLELEID,INFO/AF_ESP,INFO/AF_EXAC,INFO/AF_TGP \
  -Oz -o "$OUT" "$TMPDIR/csq.bcf"
"$TOOLS/bcftools" index -t "$OUT"
echo "records: $("$TOOLS/bcftools" view -H "$OUT" | wc -l)"
echo "with BCSQ: $("$TOOLS/bcftools" view -H -i 'INFO/BCSQ!="."' "$OUT" | wc -l)"
echo "in ClinVar: $("$TOOLS/bcftools" view -H -i 'INFO/ALLELEID!="."' "$OUT" | wc -l)"
"$TOOLS/bcftools" --version | head -1
