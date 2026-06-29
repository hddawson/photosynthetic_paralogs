#!/usr/bin/env bash
#
# find_duplicates.sh
#
# Pipeline to find duplicated genes in one annotated plant genome,
# ignoring the *origin* of the duplication (WGD / tandem / etc. are all
# lumped together). The logic is:
#
#   genome FASTA + GFF
#     -> extract one representative protein per gene  (gffread + python)
#     -> all-vs-all protein similarity search         (DIAMOND)
#     -> group genes into families (connected components of the hit graph)
#     -> per-gene table: family_id, family_size, is_duplicated
#
# A gene is called "duplicated" if its family contains >1 gene.
# family_size (count of paralogs + itself) is the response variable you will
# later regress against gene categories (cTP? GO term? photosynthetic?).
#
# Usage:
#   ./find_duplicates.sh <genome.fa[.gz]> <annotation.gff> <output_dir> [threads]
#
# This intentionally does NOT classify duplication type. The DupGen_finder
# step will be wired in later for that.

set -euo pipefail

# ---- tool paths on this machine (edit here if they move) --------------------
gffread_executable="/programs/gffread-0.9.12/gffread/gffread"
diamond_executable="/programs/diamond/diamond"
# the two helper scripts are expected to sit next to this script
script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prepare_proteins_script="${script_directory}/prepare_representative_proteins.py"
cluster_script="${script_directory}/cluster_into_families.py"
harmonize_script="${script_directory}/harmonize_gff_seqids.py"

# ---- arguments --------------------------------------------------------------
genome_fasta_input="${1:?need genome FASTA (.fa or .fa.gz)}"
annotation_gff="${2:?need annotation GFF}"
output_directory="${3:?need output directory}"
number_of_threads="${4:-8}"

# fail early and loudly if inputs or tools are missing
for required_file in "$genome_fasta_input" "$annotation_gff" \
                     "$prepare_proteins_script" "$cluster_script" \
                     "$harmonize_script" \
                     "$gffread_executable" "$diamond_executable"; do
  [ -e "$required_file" ] || { echo "ERROR: missing: $required_file" >&2; exit 1; }
done

mkdir -p "$output_directory"

# ---- step 1: make sure the genome FASTA is uncompressed ---------------------
# gffread cannot read a gzipped FASTA, so decompress into the output dir if needed.
if [[ "$genome_fasta_input" == *.gz ]]; then
  genome_fasta="${output_directory}/genome.fa"
  echo "[1/5] decompressing genome -> $genome_fasta"
  gunzip -c "$genome_fasta_input" > "$genome_fasta"
else
  genome_fasta="$genome_fasta_input"
  echo "[1/5] genome already uncompressed"
fi

# ---- step 1b: make the GFF's seqids match the FASTA's -----------------------
# GFFs name chromosomes by accession (OZ408683.1) but FASTAs use Chr1.. ; gffread
# needs them to agree. This infers the map by sequence length and aborts if the
# GFF turns out to be a different assembly than the FASTA.
harmonized_gff="${output_directory}/annotation.seqid_matched.gff"
echo "[1b/5] matching GFF seqids to FASTA"
python3 "$harmonize_script" \
  --gff "$annotation_gff" \
  --fasta "$genome_fasta" \
  --out-gff "$harmonized_gff"

# ---- step 2: extract every protein (one per transcript) ---------------------
# gffread -y writes a protein sequence per CODING transcript, headed by the
# transcript ID. Isoforms are collapsed in step 3, not here.
all_proteins_fasta="${output_directory}/all_transcript_proteins.faa"
echo "[2/5] extracting proteins with gffread"
"$gffread_executable" -y "$all_proteins_fasta" -g "$genome_fasta" "$harmonized_gff"
[ -s "$all_proteins_fasta" ] || { echo "ERROR: gffread produced no proteins (does the GFF have CDS features?)" >&2; exit 1; }

# ---- step 3: keep the longest protein per gene ------------------------------
# Collapses isoforms so that a gene's own splice variants don't masquerade as
# paralogs. Output headers are GENE ids (not transcript ids).
representative_proteins_fasta="${output_directory}/representative_proteins.faa"
transcript_to_gene_map="${output_directory}/transcript_to_gene.tsv"
echo "[3/5] collapsing isoforms to longest protein per gene"
python3 "$prepare_proteins_script" \
  --gff "$harmonized_gff" \
  --proteins "$all_proteins_fasta" \
  --out-fasta "$representative_proteins_fasta" \
  --out-map "$transcript_to_gene_map"

# ---- step 4: all-vs-all protein search with DIAMOND -------------------------
# Build a database from the representative proteins, then search them against
# themselves. We request query/subject coverage columns so the clustering step
# can require real alignments (not just a shared domain) before linking genes.
diamond_database="${output_directory}/representative_proteins.dmnd"
all_vs_all_hits="${output_directory}/all_vs_all_hits.tsv"
echo "[4/5] DIAMOND all-vs-all"
"$diamond_executable" makedb \
  --in "$representative_proteins_fasta" \
  --db "$diamond_database" \
  --threads "$number_of_threads" --quiet
"$diamond_executable" blastp \
  --query "$representative_proteins_fasta" \
  --db "$diamond_database" \
  --out "$all_vs_all_hits" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
  --evalue 1e-5 \
  --max-target-seqs 500 \
  --threads "$number_of_threads" --quiet

# ---- step 5: cluster into families and write the per-gene table -------------
per_gene_table="${output_directory}/per_gene_duplication.tsv"
families_table="${output_directory}/gene_families.tsv"
echo "[5/5] clustering into gene families"
python3 "$cluster_script" \
  --hits "$all_vs_all_hits" \
  --gene-fasta "$representative_proteins_fasta" \
  --out-per-gene "$per_gene_table" \
  --out-families "$families_table"

echo "DONE."
echo "  per-gene duplication table : $per_gene_table"
echo "  gene family membership     : $families_table"