#!/usr/bin/env python3
"""
Build a gene -> category table from the reference protein FASTAs fed into the
gene-scanning pipeline.

Categories are assigned by filename:
  - the housekeeping file              -> "null"
  - each PDB structure reference file  -> that structure's PDB ID (e.g. "9LK4")

The output gene_name matches what appears in hits_summary.tsv: the first
whitespace token of the FASTA header, with "|" and spaces replaced by "_"
(the same sanitisation the scan applies), so the table joins cleanly.

OUTPUT:
    data/reference_protein_sequences/gene_categories.tsv
        columns: gene_name, category, source_file, original_id
"""

import glob
import os
import re
import csv

# ----------------------------- configuration ------------------------------- #

reference_proteins_directory = "data/reference_protein_sequences"
housekeeping_fasta_path      = "data/housekeeping_candidates_Zm.peptides.faa"
output_tsv_path              = "data/reference_protein_sequences/gene_categories.tsv"

# --------------------------------------------------------------------------- #


def sanitise_gene_name(raw_header_id):
    """
    Reproduce the scan's gene-name sanitisation so categories join to hits:
    take the first whitespace token, then replace "|" and spaces with "_".
    e.g. "sp|P83755|PSBA_ARATH" -> "sp_P83755_PSBA_ARATH"
    """
    first_token = raw_header_id.split()[0]
    return first_token.replace("|", "_").replace(" ", "_")


def category_from_filename(fasta_path):
    """
    Decide a category label from the file name.
    The housekeeping file is "null"; a PDB reference file
    "<PDBID>_reference_proteins.faa" is labelled with its PDB ID.
    """
    filename = os.path.basename(fasta_path)

    # Housekeeping file -> null category.
    if os.path.abspath(fasta_path) == os.path.abspath(housekeeping_fasta_path):
        return "null"

    # PDB reference file: take the leading token before "_reference_proteins".
    pdb_match = re.match(r"(.+?)_reference_proteins\.faa$", filename)
    if pdb_match:
        return pdb_match.group(1)

    # Fall back to the bare filename (without extension) so nothing is silently
    # miscategorised; this is visible in the output for checking.
    return os.path.splitext(filename)[0]


def read_gene_ids_from_fasta(fasta_path):
    """Yield (sanitised_gene_name, original_header_id) for each sequence."""
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                original_id = line[1:].strip()
                yield sanitise_gene_name(original_id), original_id


def main():
    # Collect every reference FASTA the scan would use.
    fasta_extensions = ("*.fa", "*.faa", "*.fasta")
    reference_fasta_paths = []
    for extension in fasta_extensions:
        reference_fasta_paths.extend(
            glob.glob(os.path.join(reference_proteins_directory, extension))
        )
    # Include the housekeeping file even if it lives outside the reference dir.
    if os.path.exists(housekeeping_fasta_path) \
            and housekeeping_fasta_path not in reference_fasta_paths:
        reference_fasta_paths.append(housekeeping_fasta_path)

    reference_fasta_paths = sorted(set(reference_fasta_paths))
    assert len(reference_fasta_paths) > 0, \
        f"No reference FASTA files found in {reference_proteins_directory}"

    category_rows = []
    seen_gene_names = set()

    for fasta_path in reference_fasta_paths:
        category = category_from_filename(fasta_path)
        for gene_name, original_id in read_gene_ids_from_fasta(fasta_path):
            # Flag (don't silently drop) any gene appearing in two files.
            assert gene_name not in seen_gene_names, (
                f"Gene {gene_name} appears in more than one reference file; "
                "categories would be ambiguous."
            )
            seen_gene_names.add(gene_name)
            category_rows.append({
                "gene_name":   gene_name,
                "category":    category,
                "source_file": os.path.basename(fasta_path),
                "original_id": original_id,
            })

    assert len(category_rows) > 0, "No sequences read from any reference FASTA."

    os.makedirs(os.path.dirname(output_tsv_path), exist_ok=True)
    with open(output_tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(
            tsv_file,
            fieldnames=["gene_name", "category", "source_file", "original_id"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(category_rows)

    # Console summary: how many genes per category.
    counts_by_category = {}
    for row in category_rows:
        counts_by_category[row["category"]] = \
            counts_by_category.get(row["category"], 0) + 1

    print(f"Wrote {len(category_rows)} genes to {output_tsv_path}")
    print("Genes per category:")
    for category, count in sorted(counts_by_category.items()):
        print(f"  {category}: {count}")


if __name__ == "__main__":
    main()
