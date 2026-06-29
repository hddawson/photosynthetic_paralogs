#!/usr/bin/env python3
"""
Bundle all per-genome rbcS hit FASTAs into a single FASTA file and
copy it to ~/transit/ for transfer to the GPU server.

Run this after find_rbcS_copies.py has finished.

The per-genome files (rbcS_results/<accession>_rbcS_hits.fna) are simply
concatenated — Tiberius takes one multi-sequence FASTA and annotates each
sequence independently.

OUTPUT:
    rbcS_results/all_rbcS_hits_bundled.fna   (local copy)
    ~/transit/all_rbcS_hits_bundled.fna      (copy for transfer)
"""

import glob
import os
import shutil

# ----------------------------- configuration ------------------------------- #

# Directory written by find_rbcS_copies.py.
results_directory = "rbcS_results"

# Name of the bundled output file.
bundled_fasta_filename = "all_rbcS_hits_bundled.fna"

# Transit directory on this machine for cross-server file transfer.
transit_directory = os.path.expanduser("~/transit")

# --------------------------------------------------------------------------- #


def count_sequences_in_fasta(fasta_path):
    """Return the number of sequences (header lines) in a FASTA file."""
    count = 0
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                count += 1
    return count


def main():
    per_genome_fasta_paths = sorted(
        glob.glob(os.path.join(results_directory, "*_rbcS_hits.fna"))
    )
    assert len(per_genome_fasta_paths) > 0, (
        f"No *_rbcS_hits.fna files found in {results_directory}. "
        "Run find_rbcS_copies.py first."
    )
    print(f"Found {len(per_genome_fasta_paths)} per-genome FASTA files.")

    bundled_fasta_path = os.path.join(results_directory, bundled_fasta_filename)

    total_sequences = 0
    with open(bundled_fasta_path, "w") as bundled_fasta:
        for fasta_path in per_genome_fasta_paths:
            sequence_count = count_sequences_in_fasta(fasta_path)
            total_sequences += sequence_count
            print(f"  {os.path.basename(fasta_path):50s} {sequence_count} sequences")
            with open(fasta_path) as individual_fasta:
                # Write content, ensuring a blank line between files
                # doesn't accidentally merge the last and first sequences.
                content = individual_fasta.read()
                if not content.endswith("\n"):
                    content += "\n"
                bundled_fasta.write(content)

    assert count_sequences_in_fasta(bundled_fasta_path) == total_sequences, \
        "Sequence count in bundled file doesn't match sum of inputs."

    print(f"\nBundled {total_sequences} sequences -> {bundled_fasta_path}")

    # Copy to transit directory.
    os.makedirs(transit_directory, exist_ok=True)
    transit_destination = os.path.join(transit_directory, bundled_fasta_filename)
    shutil.copy2(bundled_fasta_path, transit_destination)
    print(f"Copied to: {transit_destination}")
    print(
        f"\nOn the GPU server, run:\n"
        f"  python tiberius.py --genome ~/transit/{bundled_fasta_filename} \\\n"
        f"    --model_cfg angiosperms \\\n"
        f"    --out rbcS_hits_tiberius.gtf \\\n"
        f"    --codingseq rbcS_hits_CDS.fna \\\n"
        f"    --protseq rbcS_hits_peps.faa\n"
        f"\nThen copy rbcS_hits_tiberius.gtf, rbcS_hits_CDS.fna, rbcS_hits_peps.faa "
        f"back to ~/transit/"
    )


if __name__ == "__main__":
    main()
