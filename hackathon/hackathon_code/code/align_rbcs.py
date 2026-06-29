#!/usr/bin/env python3
"""
5_align_passing_sequences.py

Extract amino acid sequences that pass both QC filters from the QC table,
write them to a FASTA, and align with MAFFT.

USAGE:
    python 5_align_passing_sequences.py

OUTPUT:
    results/rbcs_qc/rbcs_passing.faa          — unaligned FASTA of passing seqs
    results/rbcs_qc/rbcs_passing_aligned.faa  — MAFFT multiple sequence alignment
"""

import csv
import subprocess
from pathlib import Path


# ----------------------------- configuration ------------------------------- #

qc_table_path          = "results/rbcs_qc/rbcs_qc_table.tsv"
output_directory       = "results/rbcs_qc"
unaligned_fasta_path   = "results/rbcs_qc/rbcs_passing.faa"
aligned_fasta_path     = "results/rbcs_qc/rbcs_passing_aligned.faa"

mafft_binary           = "/programs/mafft/bin/mafft"

# Number of CPU threads for MAFFT.
mafft_threads          = 8

# --------------------------------------------------------------------------- #


def main():
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    # --- extract passing sequences from QC table ---
    passing_sequences = []
    with open(qc_table_path) as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            aln_passed  = row["aln_filter_passed"]  == "True"
            fold_passed = row["fold_filter_passed"] == "True"
            if aln_passed and fold_passed:
                passing_sequences.append({
                    "gene_id":          row["gene_id"],
                    "amino_acid_seq":   row["amino_acid_seq"],
                })

    assert len(passing_sequences) > 0, \
        "No sequences passed both filters — check your QC table."
    print(f"[1] sequences passing both filters: {len(passing_sequences)}")

    # --- write unaligned FASTA ---
    with open(unaligned_fasta_path, "w") as fasta_file:
        for entry in passing_sequences:
            fasta_file.write(f">{entry['gene_id']}\n")
            sequence = entry["amino_acid_seq"]
            for i in range(0, len(sequence), 80):
                fasta_file.write(sequence[i:i+80] + "\n")
    print(f"[2] unaligned FASTA written to {unaligned_fasta_path}")

    # --- run MAFFT ---
    # --auto   : MAFFT picks the best strategy for the input size
    # --thread : number of CPU threads
    # output is written to stdout and redirected to the aligned FASTA
    print(f"[3] running MAFFT...")
    with open(aligned_fasta_path, "w") as aligned_file:
        completed = subprocess.run(
            [mafft_binary, "--auto", "--thread", str(mafft_threads),
             unaligned_fasta_path],
            stdout=aligned_file,
            stderr=subprocess.PIPE,
            text=True,
        )

    assert completed.returncode == 0, \
        f"MAFFT failed:\n{completed.stderr}"

    print(f"[done] alignment written to {aligned_fasta_path}")


if __name__ == "__main__":
    main()
