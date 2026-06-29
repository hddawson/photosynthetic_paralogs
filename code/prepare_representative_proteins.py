#!/usr/bin/env python3
"""
prepare_representative_proteins.py

Collapse isoforms: keep one protein per GENE (the longest), so that a gene's
own splice variants are not later mistaken for paralogs.

Inputs:
  --gff        annotation GFF (used only to map transcript -> gene)
  --proteins   gffread -y output (one protein per transcript, headed by tx id)
Outputs:
  --out-fasta  representative proteins, headers are GENE ids
  --out-map    tab-separated transcript_id <TAB> gene_id (for your reference)

Standard library only.
"""

import argparse


def parse_gff_attributes(attributes_field):
    """Turn a GFF column-9 string 'ID=x;Parent=y;...' into a dict {ID:x,Parent:y}."""
    attribute_dictionary = {}
    for key_value_pair in attributes_field.strip().strip(";").split(";"):
        if "=" in key_value_pair:
            key, value = key_value_pair.split("=", 1)
            attribute_dictionary[key.strip()] = value.strip()
    return attribute_dictionary


def build_transcript_to_gene_map(gff_path):
    """
    Read transcript (mRNA) lines from the GFF and map each transcript ID to its
    parent gene ID. A transcript's column-9 carries ID=<transcript> and
    Parent=<gene>.
    """
    transcript_to_gene = {}
    transcript_feature_types = {"mRNA", "transcript"}
    with open(gff_path) as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue
            feature_type = columns[2]
            if feature_type not in transcript_feature_types:
                continue
            attributes = parse_gff_attributes(columns[8])
            transcript_id = attributes.get("ID")
            gene_id = attributes.get("Parent")
            if transcript_id and gene_id:
                transcript_to_gene[transcript_id] = gene_id
    # if this is empty the GFF doesn't use the ID/Parent convention we assume
    assert len(transcript_to_gene) > 0, (
        "no transcript->gene links found; check that the GFF has mRNA lines "
        "with ID= and Parent= attributes"
    )
    return transcript_to_gene


def read_fasta(fasta_path):
    """Yield (header_first_token, sequence) pairs from a FASTA file."""
    current_header = None
    current_sequence_chunks = []
    with open(fasta_path) as fasta_handle:
        for line in fasta_handle:
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(current_sequence_chunks)
                # header id is the first whitespace-delimited token after '>'
                current_header = line[1:].split()[0]
                current_sequence_chunks = []
            else:
                current_sequence_chunks.append(line.strip())
    if current_header is not None:
        yield current_header, "".join(current_sequence_chunks)


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--gff", required=True)
    argument_parser.add_argument("--proteins", required=True)
    argument_parser.add_argument("--out-fasta", required=True)
    argument_parser.add_argument("--out-map", required=True)
    arguments = argument_parser.parse_args()

    transcript_to_gene = build_transcript_to_gene_map(arguments.gff)

    # For each gene, remember the longest protein seen across its transcripts.
    longest_protein_sequence_for_gene = {}
    number_of_transcripts_without_gene = 0

    for transcript_id, protein_sequence in read_fasta(arguments.proteins):
        # if a protein header isn't in the GFF map, treat that protein as its
        # own gene rather than dropping it (and count it so we can warn)
        gene_id = transcript_to_gene.get(transcript_id)
        if gene_id is None:
            gene_id = transcript_id
            number_of_transcripts_without_gene += 1
        previous_best = longest_protein_sequence_for_gene.get(gene_id)
        if previous_best is None or len(protein_sequence) > len(previous_best):
            longest_protein_sequence_for_gene[gene_id] = protein_sequence

    assert len(longest_protein_sequence_for_gene) > 0, "no proteins read"

    # write one representative protein per gene, headed by the gene id
    with open(arguments.out_fasta, "w") as fasta_output:
        for gene_id, protein_sequence in longest_protein_sequence_for_gene.items():
            fasta_output.write(f">{gene_id}\n{protein_sequence}\n")

    # write the transcript -> gene mapping for the record
    with open(arguments.out_map, "w") as map_output:
        map_output.write("transcript_id\tgene_id\n")
        for transcript_id, gene_id in transcript_to_gene.items():
            map_output.write(f"{transcript_id}\t{gene_id}\n")

    print(f"  genes with a representative protein : {len(longest_protein_sequence_for_gene)}")
    if number_of_transcripts_without_gene:
        print(f"  WARNING: {number_of_transcripts_without_gene} proteins had no "
              f"GFF gene parent and were treated as their own gene")


if __name__ == "__main__":
    main()
