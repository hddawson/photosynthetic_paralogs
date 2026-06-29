#!/usr/bin/env bash
#
# run_all_genomes.sh
#
# Run find_duplicates.sh on every genome that pair_genomes.sh resolved to a real
# FASTA. Each genome runs independently: if one fails (e.g. an assembly mismatch
# caught at step 1b), it is logged and the loop continues to the next.
#
# Usage:
#   ./run_all_genomes.sh <pairs.tsv> <gff_dir> <results_dir> [threads]
#
# pairs.tsv is the output of pair_genomes.sh (columns: gff_file, genome_path, status).

set -euo pipefail

pairs_table="${1:?need pairs.tsv from pair_genomes.sh}"
gff_directory="${2:?need gff directory}"
results_directory="${3:?need results directory}"
number_of_threads="${4:-8}"

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
find_duplicates_script="${script_directory}/find_duplicates.sh"

mkdir -p "$results_directory"
run_log="${results_directory}/run_status.tsv"
printf "sample\tstatus\n" > "$run_log"

# iterate only over rows pair_genomes.sh marked OK (skip the header and MISSING)
while IFS=$'\t' read -r gff_file genome_path status; do
  [ "$status" = "OK" ] || continue
  sample=${gff_file%_GeneCAD_final.gff}
  echo "=================== $sample ==================="

  # run one genome; capture success/failure without aborting the whole batch
  if bash "$find_duplicates_script" \
        "$genome_path" \
        "${gff_directory}/${gff_file}" \
        "${results_directory}/${sample}" \
        "$number_of_threads" \
        > "${results_directory}/${sample}.log" 2>&1; then
    printf "%s\tOK\n" "$sample" | tee -a "$run_log"
  else
    printf "%s\tFAILED\n" "$sample" | tee -a "$run_log"
    echo "  see ${results_directory}/${sample}.log for the error"
  fi
done < "$pairs_table"

echo ""
echo "batch finished. summary:"
# count outcomes (skip header)
awk -F'\t' 'NR>1{c[$2]++} END{for(s in c) print "  "s": "c[s]}' "$run_log"
echo "full status table: $run_log"
