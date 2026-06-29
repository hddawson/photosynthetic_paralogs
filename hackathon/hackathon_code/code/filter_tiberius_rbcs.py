#!/usr/bin/env python3
"""
Parse Tiberius peptide output and write a TSV of rbcS candidates.

Filtering: length only. The EPYCIK motif filter was too strict — real rbcS
proteins from diverse angiosperms do not all carry that exact variant.
Sequences shorter than minimum_peptide_length_aa are discarded as obvious
fragments; everything else goes into the output TSV for inspection in R.

Input:  ~/transit/rbcS_hits_peps.faa   (Tiberius --protseq output)
Output: rbcS_results/rbcS_peptides.tsv
            columns: peptide_id, assembly_accession, species, length_aa, sequence
"""

import os
import csv

# ----------------------------- configuration ------------------------------- #

transit_directory        = os.path.expanduser("~/transit")
tiberius_peptide_fasta   = os.path.join(transit_directory, "rbcS_hits_peps.faa")
results_directory        = "rbcS_results"

# Minimum peptide length to keep. Removes only obvious truncation artifacts.
# The shortest real rbcS mature domain is ~100 aa; transit peptide adds ~50 aa.
minimum_peptide_length_aa = 80

# --------------------------------------------------------------------------- #


def parse_peptide_id(header_line):
    """
    Extract assembly accession and species from a Tiberius peptide header.

    Tiberius header format (from our bundled input):
        g1|g10.t1|GCA_018105755.1_rbcS_copy6:1479-2328(-)

    The GCA accession is embedded in the third pipe-delimited field, before
    the underscore that precedes '_rbcS_copy'.
    """
    # Full ID is everything after '>'
    peptide_id = header_line.lstrip(">").split()[0]

    # Pull the GCA accession from the third |-delimited token.
    fields = peptide_id.split("|")
    assembly_accession = "unknown"
    if len(fields) >= 3:
        # e.g. "GCA_018105755.1_rbcS_copy6:1479-2328(-)"
        third_field = fields[2]
        # GCA accession ends before the first underscore after the version dot
        # pattern: GCA_XXXXXXXXX.V_rbcS...
        parts = third_field.split("_rbcS")
        if parts:
            assembly_accession = parts[0]

    return peptide_id, assembly_accession


def read_fasta_records(fasta_path):
    """
    Yield (header_line, sequence_string) tuples from a FASTA file.
    header_line includes the leading '>'.
    """
    current_header = None
    current_parts  = []
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(current_parts)
                current_header = line
                current_parts  = []
            elif line:
                current_parts.append(line)
    if current_header is not None:
        yield current_header, "".join(current_parts)


def main():
    assert os.path.exists(tiberius_peptide_fasta), (
        f"Tiberius peptide file not found: {tiberius_peptide_fasta}\n"
        "Copy rbcS_hits_peps.faa from the GPU server to ~/transit/ first."
    )

    os.makedirs(results_directory, exist_ok=True)
    output_tsv_path = os.path.join(results_directory, "rbcS_peptides.tsv")

    total_read     = 0
    total_kept     = 0
    total_too_short = 0

    with open(output_tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(
            tsv_file,
            fieldnames=["peptide_id", "assembly_accession", "length_aa", "sequence"],
            delimiter="\t",
        )
        writer.writeheader()

        for header_line, sequence in read_fasta_records(tiberius_peptide_fasta):
            total_read += 1
            length_aa   = len(sequence)

            if length_aa < minimum_peptide_length_aa:
                total_too_short += 1
                continue

            peptide_id, assembly_accession = parse_peptide_id(header_line)

            writer.writerow({
                "peptide_id":         peptide_id,
                "assembly_accession": assembly_accession,
                "length_aa":          length_aa,
                "sequence":           sequence,
            })
            total_kept += 1

    print(f"Read:       {total_read} peptides")
    print(f"Too short:  {total_too_short}  (< {minimum_peptide_length_aa} aa)")
    print(f"Written:    {total_kept} peptides -> {output_tsv_path}")


if __name__ == "__main__":
    main()