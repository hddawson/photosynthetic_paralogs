#!/usr/bin/env python3
"""
extract_gene_positions.py

Write one row per gene: gene_id, chromosome, start, end. Positions come from the
gene's transcript lines (gene_id = mRNA Parent, span = min start..max end), so
gene ids match representative_proteins.faa and the category labels exactly.

Inputs:
  --gff          harmonized annotation GFF (seqids match the FASTA)
  --genes-from   representative_proteins.faa (defines which genes to include)
Output:
  --out          gene_id <TAB> chromosome <TAB> start <TAB> end
"""

import argparse


def read_gene_ids(fasta_path):
    gene_ids = []
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                gene_ids.append(line[1:].split()[0])
    assert gene_ids, f"no gene ids in {fasta_path}"
    return gene_ids


def parse_gff_attributes(attributes_field):
    attribute_dictionary = {}
    for key_value_pair in attributes_field.strip().strip(";").split(";"):
        if "=" in key_value_pair:
            key, value = key_value_pair.split("=", 1)
            attribute_dictionary[key.strip()] = value.strip()
    return attribute_dictionary


def read_gene_positions(gff_path):
    chromosome_of_gene, start_of_gene, end_of_gene = {}, {}, {}
    transcript_feature_types = {"mRNA", "transcript"}
    with open(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9 or columns[2] not in transcript_feature_types:
                continue
            gene_id = parse_gff_attributes(columns[8]).get("Parent")
            if gene_id is None:
                continue
            chromosome = columns[0]
            start_coordinate = int(columns[3])
            end_coordinate = int(columns[4])
            chromosome_of_gene.setdefault(gene_id, chromosome)
            if gene_id not in start_of_gene or start_coordinate < start_of_gene[gene_id]:
                start_of_gene[gene_id] = start_coordinate
            if gene_id not in end_of_gene or end_coordinate > end_of_gene[gene_id]:
                end_of_gene[gene_id] = end_coordinate
    assert chromosome_of_gene, f"no transcript positions found in {gff_path}"
    return chromosome_of_gene, start_of_gene, end_of_gene


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--genes-from", required=True)
    parser.add_argument("--out", required=True)
    arguments = parser.parse_args()

    genes_to_include = read_gene_ids(arguments.genes_from)
    chromosome_of_gene, start_of_gene, end_of_gene = read_gene_positions(arguments.gff)

    number_written, number_missing = 0, 0
    with open(arguments.out, "w") as out_handle:
        out_handle.write("gene_id\tchromosome\tstart\tend\n")
        for gene_id in genes_to_include:
            if gene_id not in chromosome_of_gene:
                number_missing += 1
                continue
            out_handle.write(f"{gene_id}\t{chromosome_of_gene[gene_id]}\t"
                             f"{start_of_gene[gene_id]}\t{end_of_gene[gene_id]}\n")
            number_written += 1

    print(f"  wrote {number_written} gene positions")
    if number_missing:
        print(f"  WARNING: {number_missing} proteins had no GFF position")


if __name__ == "__main__":
    main()
