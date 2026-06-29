#!/usr/bin/env bash
#
# label_gene_categories.sh
#
# Assign each gene of one genome a category (photosynthetic complex or
# housekeeping) by DIAMOND-searching its representative proteins against the
# reference protein set, then taking the best hit's category.
# Writes results/<sample>/gene_category.tsv  (gene_id, category, ...).
#
# Usage:
#   ./label_gene_categories.sh <sample> [results_root] [threads]

set -euo pipefail

sample_name="${1:?need sample name}"
results_root="${2:-results}"
number_of_threads="${3:-8}"

diamond_executable="/programs/diamond/diamond"
script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
assign_category_script="${script_directory}/assign_best_category.py"

# user-provided reference material
reference_directory="data/reference_protein_sequences"
housekeeping_fasta="data/housekeeping_candidates_Zm.peptides.faa"
gene_categories_table="${reference_directory}/gene_categories.tsv"

target_proteins="${results_root}/${sample_name}/representative_proteins.faa"

for required_file in "$target_proteins" "$assign_category_script" \
                     "$gene_categories_table" "$diamond_executable"; do
  [ -e "$required_file" ] || { echo "ERROR: missing: $required_file" >&2; exit 1; }
done

# build the combined reference DB once, cached and reused for every species
reference_shared_directory="${results_root}/_photosynthesis_reference"
mkdir -p "$reference_shared_directory"
combined_reference_fasta="${reference_shared_directory}/all_reference_proteins.faa"
combined_reference_database="${reference_shared_directory}/all_reference_proteins"

if [ ! -s "${combined_reference_database}.dmnd" ]; then
  echo "building combined reference DB (once)"
  # concatenate every reference FASTA plus the housekeeping file
  : > "$combined_reference_fasta"
  for reference_fasta in "$reference_directory"/*.fa "$reference_directory"/*.faa \
                         "$reference_directory"/*.fasta "$housekeeping_fasta"; do
    [ -e "$reference_fasta" ] && cat "$reference_fasta" >> "$combined_reference_fasta"
  done
  [ -s "$combined_reference_fasta" ] || { echo "ERROR: no reference FASTAs found" >&2; exit 1; }
  "$diamond_executable" makedb --in "$combined_reference_fasta" \
    --db "$combined_reference_database" --threads "$number_of_threads" --quiet
fi

work_directory="${results_root}/${sample_name}/category"
mkdir -p "$work_directory"
hits_file="${work_directory}/reference_hits.tsv"

echo "DIAMOND ${sample_name} proteins vs reference set"
"$diamond_executable" blastp -q "$target_proteins" -d "$combined_reference_database" \
  -o "$hits_file" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
  --evalue 1e-8 --max-target-seqs 5 --threads "$number_of_threads" --quiet

python3 "$assign_category_script" \
  --hits "$hits_file" --categories "$gene_categories_table" \
  --out "${results_root}/${sample_name}/gene_category.tsv"

echo "DONE: ${results_root}/${sample_name}/gene_category.tsv"
