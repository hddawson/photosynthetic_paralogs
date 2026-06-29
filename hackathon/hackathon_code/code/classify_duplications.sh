#!/usr/bin/env bash
#
# classify_duplications.sh
#
# Run DupGen_finder for one target species against a fixed outgroup, reusing the
# representative proteins and harmonized GFF that find_duplicates.sh produced.
# Adds results/<target>/per_gene_duplication_type.tsv  (gene_id, duplication_type)
# with types: wgd / tandem / proximal / transposed / dispersed.
#
# Usage:
#   ./classify_duplications.sh <target_sample> <outgroup_sample> [results_root] [threads]

set -euo pipefail

target_sample="${1:?need target sample name}"
outgroup_sample="${2:?need outgroup sample name}"
results_root="${3:-results}"
number_of_threads="${4:-8}"

# ---- tool paths (edit if they move) -----------------------------------------
diamond_executable="/programs/diamond/diamond"
dupgen_finder_executable="${HOME}/dupgen/DupGen_finder/DupGen_finder.pl"
script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
make_position_file_script="${script_directory}/make_mcscanx_gff.py"
summarize_types_script="${script_directory}/summarize_dupgen_types.py"

# ---- inputs produced by find_duplicates.sh ----------------------------------
target_proteins="${results_root}/${target_sample}/representative_proteins.faa"
target_gff="${results_root}/${target_sample}/annotation.seqid_matched.gff"
outgroup_proteins="${results_root}/${outgroup_sample}/representative_proteins.faa"
outgroup_gff="${results_root}/${outgroup_sample}/annotation.seqid_matched.gff"

for required_file in "$target_proteins" "$target_gff" \
                     "$outgroup_proteins" "$outgroup_gff" \
                     "$make_position_file_script" "$summarize_types_script" \
                     "$diamond_executable" "$dupgen_finder_executable"; do
  [ -e "$required_file" ] || { echo "ERROR: missing: $required_file" >&2; exit 1; }
done

work_directory="${results_root}/${target_sample}/dupgen"
data_directory="${work_directory}/data"
database_directory="${work_directory}/db"
output_directory="${work_directory}/out"
mkdir -p "$data_directory" "$database_directory" "$output_directory"

# the outgroup's prefixed proteins and position file are identical for every
# target, so build them once and cache under the outgroup's own folder
outgroup_shared_directory="${results_root}/${outgroup_sample}/dupgen_shared"
mkdir -p "$outgroup_shared_directory"
outgroup_position_file="${outgroup_shared_directory}/outgroup.gff"
outgroup_prefixed_proteins="${outgroup_shared_directory}/outgroup.prefixed.faa"
outgroup_database="${outgroup_shared_directory}/outgroup"

echo "[1/5] target gene positions"
python3 "$make_position_file_script" \
  --gff "$target_gff" --genes-from "$target_proteins" \
  --chrom-tag ta --out "${data_directory}/${target_sample}.gff"

echo "[2/5] outgroup gene positions + prefixed protein DB (cached)"
if [ ! -s "$outgroup_position_file" ]; then
  python3 "$make_position_file_script" \
    --gff "$outgroup_gff" --genes-from "$outgroup_proteins" \
    --chrom-tag og --gene-prefix OG_ --out "$outgroup_position_file"
fi
if [ ! -s "${outgroup_database}.dmnd" ]; then
  # prefix every outgroup protein header with OG_ so BLAST subject ids match
  # the prefixed gene ids in the outgroup position file
  awk '/^>/{print ">OG_"substr($0,2); next}{print}' \
    "$outgroup_proteins" > "$outgroup_prefixed_proteins"
  "$diamond_executable" makedb --in "$outgroup_prefixed_proteins" \
    --db "$outgroup_database" --threads "$number_of_threads" --quiet
fi

# combined position file = target + outgroup (distinct tags, no id collisions)
cat "${data_directory}/${target_sample}.gff" "$outgroup_position_file" \
  > "${data_directory}/${target_sample}_${outgroup_sample}.gff"

echo "[3/5] DIAMOND self (target vs target)"
# 12-column outfmt 6 and top-5 hits at e<=1e-10 is the DupGen_finder convention
"$diamond_executable" makedb --in "$target_proteins" \
  --db "${database_directory}/target" --threads "$number_of_threads" --quiet
"$diamond_executable" blastp -q "$target_proteins" -d "${database_directory}/target" \
  -o "${data_directory}/${target_sample}.blast" \
  --outfmt 6 --evalue 1e-10 --max-target-seqs 5 \
  --threads "$number_of_threads" --quiet

echo "[4/5] DIAMOND target vs outgroup"
"$diamond_executable" blastp -q "$target_proteins" -d "$outgroup_database" \
  -o "${data_directory}/${target_sample}_${outgroup_sample}.blast" \
  --outfmt 6 --evalue 1e-10 --max-target-seqs 5 \
  --threads "$number_of_threads" --quiet

echo "[5/5] DupGen_finder + per-gene type summary"
perl "$dupgen_finder_executable" \
  -i "$data_directory" -t "$target_sample" -c "$outgroup_sample" \
  -o "$output_directory"

python3 "$summarize_types_script" \
  --dupgen-out-dir "$output_directory" --target "$target_sample" \
  --out "${results_root}/${target_sample}/per_gene_duplication_type.tsv"

echo "DONE: ${results_root}/${target_sample}/per_gene_duplication_type.tsv"
