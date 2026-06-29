#!/usr/bin/env bash
#
# run_all_labeling.sh
#
# Label gene categories for every species that find_duplicates.sh ran
# successfully. Each runs independently; failures are logged, loop continues.
#
# Usage:
#   ./run_all_labeling.sh [results_root] [threads]

set -euo pipefail

results_root="${1:-results}"
number_of_threads="${2:-8}"

echo "=== labeling gene categories for all species in $results_root ==="

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
label_script="${script_directory}/label_gene_categories.sh"
run_status_table="${results_root}/run_status.tsv"
[ -e "$run_status_table" ] || { echo "ERROR: missing $run_status_table" >&2; exit 1; }

labeling_log="${results_root}/labeling_status.tsv"
printf "sample\tstatus\n" > "$labeling_log"

while IFS=$'\t' read -r sample status; do
  [ "$status" = "OK" ] || continue
  echo "=================== $sample ==================="
  if bash "$label_script" "$sample" "$results_root" "$number_of_threads" \
        > "${results_root}/${sample}/labeling.log" 2>&1; then
    printf "%s\tOK\n" "$sample" | tee -a "$labeling_log"
  else
    printf "%s\tFAILED\n" "$sample" | tee -a "$labeling_log"
    echo "  see ${results_root}/${sample}/labeling.log"
  fi
done < "$run_status_table"

echo ""
echo "labeling finished. summary:"
awk -F'\t' 'NR>1{c[$2]++} END{for(s in c) print "  "s": "c[s]}' "$labeling_log"
