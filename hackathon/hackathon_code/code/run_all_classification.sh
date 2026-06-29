#!/usr/bin/env bash
#
# run_all_classification.sh
#
# Run classify_duplications.sh on every species that find_duplicates.sh ran
# successfully, EXCEPT the outgroup itself. Reads results/run_status.tsv to know
# which species succeeded. Each species runs independently; a failure is logged
# and the loop continues.
#
# Usage:
#   ./run_all_classification.sh <outgroup_sample> [results_root] [threads]

set -euo pipefail

outgroup_sample="${1:?need outgroup sample name, e.g. Ncolorata}"
results_root="${2:-results}"
number_of_threads="${3:-8}"

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
classify_script="${script_directory}/classify_duplications.sh"
run_status_table="${results_root}/run_status.tsv"
[ -e "$run_status_table" ] || { echo "ERROR: missing $run_status_table" >&2; exit 1; }

classification_log="${results_root}/classification_status.tsv"
printf "sample\tstatus\n" > "$classification_log"

# also confirm the outgroup itself completed find_duplicates.sh (we need its
# proteins and GFF as the reference)
awk -F'\t' -v og="$outgroup_sample" '$1==og && $2=="OK"{found=1} END{exit !found}' \
  "$run_status_table" \
  || { echo "ERROR: outgroup $outgroup_sample is not OK in $run_status_table" >&2; exit 1; }

# iterate over successful samples, skipping the outgroup
while IFS=$'\t' read -r sample status; do
  [ "$status" = "OK" ] || continue
  [ "$sample" = "$outgroup_sample" ] && continue
  echo "=================== $sample ==================="
  if bash "$classify_script" "$sample" "$outgroup_sample" \
        "$results_root" "$number_of_threads" \
        > "${results_root}/${sample}/dupgen.log" 2>&1; then
    printf "%s\tOK\n" "$sample" | tee -a "$classification_log"
  else
    printf "%s\tFAILED\n" "$sample" | tee -a "$classification_log"
    echo "  see ${results_root}/${sample}/dupgen.log"
  fi
done < "$run_status_table"

echo ""
echo "classification finished. summary:"
awk -F'\t' 'NR>1{c[$2]++} END{for(s in c) print "  "s": "c[s]}' "$classification_log"
