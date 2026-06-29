#!/usr/bin/env python3
"""
Copy the all_hits.fna output from find_rbcS_copies.py to ~/transit/
for transfer to the GPU server running Tiberius.

Run this after find_rbcS_copies.py has finished.
"""

import os
import shutil

# ----------------------------- configuration ------------------------------- #

# The single all-genes FASTA written by find_rbcS_copies.py.
all_hits_fasta_path = "results/gene_hits/all_hits.fna"

# Name to use in transit directory.
bundled_fasta_filename = "all_hits_bundled.fna"

# Transit directory for cross-server transfer.
transit_directory = os.path.expanduser("~/transit")

# --------------------------------------------------------------------------- #


def count_sequences_in_fasta(fasta_path):
    """Return the number of '>' header lines in a FASTA file."""
    count = 0
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                count += 1
    return count


def main():
    assert os.path.exists(all_hits_fasta_path), (
        f"all_hits.fna not found at {all_hits_fasta_path}. "
        "Run find_rbcS_copies.py first."
    )

    total_sequences = count_sequences_in_fasta(all_hits_fasta_path)
    assert total_sequences > 0, f"No sequences found in {all_hits_fasta_path}."
    print(f"Found {total_sequences} sequences in {all_hits_fasta_path}.")

    os.makedirs(transit_directory, exist_ok=True)
    transit_destination = os.path.join(transit_directory, bundled_fasta_filename)
    shutil.copy2(all_hits_fasta_path, transit_destination)
    print(f"Copied to: {transit_destination}")
    print(
        f"\nOn the GPU server, run:\n"
        f"  python tiberius.py --genome ~/transit/{bundled_fasta_filename} \\\n"
        f"    --model_cfg angiosperms \\\n"
        f"    --out all_hits_tiberius.gtf \\\n"
        f"    --codingseq all_hits_CDS.fna \\\n"
        f"    --protseq all_hits_peps.faa\n"
        f"\nThen copy all_hits_tiberius.gtf, all_hits_CDS.fna, all_hits_peps.faa "
        f"back to ~/transit/"
    )


if __name__ == "__main__":
    main()