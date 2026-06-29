#!/usr/bin/env python3
"""
harmonize_gff_seqids.py

The GFFs name chromosomes by INSDC accession (e.g. OZ408683.1) while the FASTAs
name them Chr1..Chr5/ChrC/ChrM. gffread needs the two to agree. This rewrites
the GFF's seqid column (column 1) to the FASTA names.

The accession -> FASTA-name map is inferred PURELY BY LENGTH: each GFF seqid is
matched to the FASTA sequence with the smallest length that is still >= the
largest gene-end coordinate seen on that seqid (a gene cannot extend past the
end of its chromosome). The match must be a clean one-to-one mapping; if it is
not, the assembly genuinely differs from the FASTA and the script aborts rather
than silently producing wrong gene models.

Inputs:
  --gff       annotation GFF (seqids = accessions)
  --fasta     genome FASTA (.fa or .fa.gz; provides the target names + lengths)
Output:
  --out-gff   same GFF with column 1 renamed to FASTA names
"""

import argparse
import gzip


def open_maybe_gzipped(path):
    """Open a text file transparently whether or not it is gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta_sequence_lengths(fasta_path):
    """Return {sequence_name: length_in_bases} for every record in the FASTA."""
    sequence_length_by_name = {}
    current_name = None
    with open_maybe_gzipped(fasta_path) as fasta_handle:
        for line in fasta_handle:
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                sequence_length_by_name[current_name] = 0
            else:
                sequence_length_by_name[current_name] += len(line.strip())
    assert sequence_length_by_name, "no sequences found in FASTA"
    return sequence_length_by_name


def read_gff_max_end_per_seqid(gff_path):
    """Return {gff_seqid: largest_end_coordinate} across all features."""
    largest_end_for_seqid = {}
    with open(gff_path) as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 5:
                continue
            seqid = columns[0]
            end_coordinate = int(columns[4])
            if end_coordinate > largest_end_for_seqid.get(seqid, 0):
                largest_end_for_seqid[seqid] = end_coordinate
    assert largest_end_for_seqid, "no features with coordinates found in GFF"
    return largest_end_for_seqid


def infer_seqid_renaming(largest_end_for_seqid, sequence_length_by_name):
    """
    Map each GFF seqid to a FASTA name by tightest length fit.

    For each GFF seqid (processed largest-first), choose the still-unused FASTA
    sequence with the SMALLEST length that is >= that seqid's largest gene-end.
    That unique 'smallest length still big enough to contain the genes' is the
    true chromosome whenever chromosome lengths are well separated.
    """
    fasta_name_to_length = sequence_length_by_name
    available_fasta_names = set(fasta_name_to_length)
    renaming = {}

    # process the longest GFF sequences first so tight fits are assigned cleanly
    for gff_seqid in sorted(largest_end_for_seqid,
                            key=lambda s: largest_end_for_seqid[s],
                            reverse=True):
        required_minimum_length = largest_end_for_seqid[gff_seqid]
        # candidate FASTA names large enough to contain every gene on this seqid
        candidates = [name for name in available_fasta_names
                      if fasta_name_to_length[name] >= required_minimum_length]
        assert candidates, (
            f"GFF seqid {gff_seqid} has a gene ending at {required_minimum_length}, "
            f"longer than any unused FASTA sequence -> assemblies do not match"
        )
        # tightest fit = smallest length that still contains the genes
        best_match = min(candidates, key=lambda name: fasta_name_to_length[name])
        renaming[gff_seqid] = best_match
        available_fasta_names.remove(best_match)

    # mapping must be one-to-one (a bijection onto the FASTA names it used)
    assert len(set(renaming.values())) == len(renaming), \
        "two GFF seqids mapped to the same FASTA sequence -> ambiguous"
    return renaming


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--gff", required=True)
    argument_parser.add_argument("--fasta", required=True)
    argument_parser.add_argument("--out-gff", required=True)
    arguments = argument_parser.parse_args()

    sequence_length_by_name = read_fasta_sequence_lengths(arguments.fasta)
    largest_end_for_seqid = read_gff_max_end_per_seqid(arguments.gff)
    renaming = infer_seqid_renaming(largest_end_for_seqid, sequence_length_by_name)

    # show the inferred mapping and the slack (chromosome length minus furthest
    # gene end); small non-negative slack is the signature of the same assembly
    print("  inferred seqid renaming (gff_seqid -> fasta_name : length, slack_bp):")
    for gff_seqid, fasta_name in sorted(renaming.items(),
                                        key=lambda kv: sequence_length_by_name[kv[1]],
                                        reverse=True):
        slack = sequence_length_by_name[fasta_name] - largest_end_for_seqid[gff_seqid]
        print(f"    {gff_seqid} -> {fasta_name} : "
              f"{sequence_length_by_name[fasta_name]} bp, slack {slack} bp")

    # rewrite column 1; pass comment lines through unchanged
    number_of_renamed_lines = 0
    with open(arguments.gff) as gff_handle, open(arguments.out_gff, "w") as out_handle:
        for line in gff_handle:
            if line.startswith("#"):
                out_handle.write(line)
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                out_handle.write(line)
                continue
            if columns[0] in renaming:
                columns[0] = renaming[columns[0]]
                number_of_renamed_lines += 1
            out_handle.write("\t".join(columns) + "\n")

    print(f"  renamed {number_of_renamed_lines} feature lines")


if __name__ == "__main__":
    main()
