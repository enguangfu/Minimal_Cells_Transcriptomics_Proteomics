#!/usr/bin/env bash
# gather_SI_excel.sh
# Copy the three Supplementary Data Excel workbooks (S1-S3) into this SI folder
# (HERE), next to the S4 QC archive, ready for upload to the journal.
# Sources live in their analysis folders; this only gathers the latest builds.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$(cd "$HERE/../.." && pwd)"

# S#  |  source (relative to project root)            |  destination filename (in HERE)
ROWS=(
  "S1|Syn1_Operon/operon.xlsx|Supplementary_Data_S1_operon.xlsx"
  "S2|Syn1_Corr_RNA_Proteins/syn1_omics.xlsx|Supplementary_Data_S2_syn1_omics.xlsx"
  "S3|Genome_Reduction/genome_reduction.xlsx|Supplementary_Data_S3_genome_reduction.xlsx"
)

fail=0
for row in "${ROWS[@]}"; do
  IFS='|' read -r label src dst <<< "$row"
  if [[ -f "$PROJECT/$src" ]]; then
    cp -f "$PROJECT/$src" "$HERE/$dst"
    printf 'OK    %-3s %s\n' "$label" "$dst"
  else
    printf 'MISS  %-3s source not found: %s\n' "$label" "$src" >&2
    fail=1
  fi
done

echo
echo "SI folder now holds:"
ls -1 "$HERE"/Supplementary_Data_* 2>/dev/null

exit $fail
