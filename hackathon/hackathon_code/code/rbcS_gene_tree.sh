#!/usr/bin/env bash
# align_and_tree_rbcS.sh
#
# Labels rbcS peptide sequences with species names, aligns with MAFFT,
# and infers a maximum-likelihood tree with RAxML-NG.
#
# REQUIRES:
#   - rbcS_results/rbcS_peptides.tsv   (from filter_tiberius_rbcS.py)
#   - rbcS_results/hits_summary.tsv    (from find_rbcS_copies.py)
#   - /programs/mafft/bin/mafft
#   - /programs/raxml-ng_v1.0.1/raxml-ng
#
# OUTPUT (all in rbcS_results/):
#   rbcS_labelled.faa          FASTA with species-named headers
#   rbcS_aligned.faa           MAFFT alignment
#   rbcS_tree.raxml.support    Best ML tree with bootstrap support

set -euo pipefail   # exit on error, unset variable, or pipe failure

RESULTS_DIR="rbcS_results"
PEPTIDE_TSV="${RESULTS_DIR}/rbcS_peptides.tsv"
LABELLED_FASTA="${RESULTS_DIR}/rbcS_labelled.faa"
ALIGNED_FASTA="${RESULTS_DIR}/rbcS_aligned.faa"
TREE_PREFIX="${RESULTS_DIR}/rbcS_tree"

MAFFT="/programs/mafft/bin/mafft"
RAXML="/programs/raxml-ng_v1.0.1/raxml-ng"

THREADS=4
BOOTSTRAP_REPLICATES=100
RAXML_MODEL="LG+G4"
RAXML_SEED=42

# --------------------------------------------------------------------------- #

echo "=== Step 1: checking inputs and tools ==="

for required_file in "${PEPTIDE_TSV}"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: required file not found: ${required_file}"
        exit 1
    fi
done

for required_binary in "${MAFFT}" "${RAXML}"; do
    if [[ ! -x "${required_binary}" ]]; then
        echo "ERROR: binary not found or not executable: ${required_binary}"
        exit 1
    fi
done

echo "  All inputs and tools found."

# --------------------------------------------------------------------------- #

echo ""
echo "=== Step 2: labelling sequences with species names ==="

# Join species from hits_summary.tsv onto peptides.tsv, then write FASTA.
# Label format: Genus_species_GCA_XXXXXXXXX.V_gN.tN
# Spaces in species names are replaced with underscores (required by RAxML).
python3 - << 'PYEOF'
import csv, os

results_dir  = "rbcS_results"
peptide_tsv  = os.path.join(results_dir, "rbcS_peptides.tsv")
output_fasta = os.path.join(results_dir, "rbcS_labelled.faa")

# Species names not yet available — label as accession + gene ID.
# Join species names in R using the NCBI metadata once the pipeline finishes.
sequences_written = 0
with open(peptide_tsv) as pep_file, open(output_fasta, "w") as fasta_out:
    for row in csv.DictReader(pep_file, delimiter="\t"):
        accession  = row["assembly_accession"]
        peptide_id = row["peptide_id"]
        sequence   = row["sequence"]

        # Use the short gene ID (e.g. g10.t1) from the middle pipe-field.
        short_gene_id = peptide_id.split("|")[1] if "|" in peptide_id else peptide_id
        # Replace characters RAxML dislikes in sequence labels.
        label = f"{accession}_{short_gene_id}".replace(":", "_").replace("(", "").replace(")", "")

        fasta_out.write(f">{label}\n{sequence}\n")
        sequences_written += 1

print(f"  Wrote {sequences_written} labelled sequences to {output_fasta}")
PYEOF

# --------------------------------------------------------------------------- #

echo ""
echo "=== Step 3: aligning with MAFFT ==="

"${MAFFT}" --auto --thread "${THREADS}" "${LABELLED_FASTA}" > "${ALIGNED_FASTA}"

SEQUENCE_COUNT=$(grep -c "^>" "${ALIGNED_FASTA}")
echo "  Aligned ${SEQUENCE_COUNT} sequences -> ${ALIGNED_FASTA}"

# --------------------------------------------------------------------------- #

echo ""
echo "=== Step 4: inferring ML tree with RAxML-NG ==="
echo "  Model: ${RAXML_MODEL}, bootstrap replicates: ${BOOTSTRAP_REPLICATES}"

# Remove any previous RAxML run with the same prefix to avoid the
# "output files already exist" error.
rm -f "${TREE_PREFIX}".raxml.*

"${RAXML}" \
    --all \
    --msa      "${ALIGNED_FASTA}" \
    --model    "${RAXML_MODEL}" \
    --bs-trees "${BOOTSTRAP_REPLICATES}" \
    --prefix   "${TREE_PREFIX}" \
    --threads  "${THREADS}" \
    --seed     "${RAXML_SEED}"

# --------------------------------------------------------------------------- #

echo ""
echo "=== Done ==="
echo "  Best ML tree with bootstrap support: ${TREE_PREFIX}.raxml.support"