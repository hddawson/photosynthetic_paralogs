#!/usr/bin/env python3
"""
make_mcscanx_gff.py

Write the gene-position file DupGen_finder/MCScanX expects: four tab-separated
columns  ->  chromosome_tag <TAB> gene_id <TAB> start <TAB> end

Two things this handles that matter for correctness:
  * chromosome labels are rewritten to <chrom_tag><index> (e.g. ta1, ta2) so the
    target and the outgroup get DISTINCT tags -- MCScanX uses the leading letters
    to tell genomes apart, and would merge them if they shared names like Chr01.
  * an optional gene-id prefix (e.g. OG_) namespaces the outgroup's gene IDs so
    they cannot collide with identically-named target genes in the combined file.

The gene set and IDs are taken from the representative-protein FASTA, so the
position file lists exactly the genes that are in the BLAST search (one per gene).

Inputs:
  --gff          harmonized annotation GFF (seqids already match the FASTA)
  --genes-from   representative_proteins.faa (defines which genes to include)
  --chrom-tag    two-letter genome tag for column 1 (e.g. 'ta' or 'og')
  --gene-prefix  string prepended to every gene id in column 2 (default none)
Output:
  --out          MCScanX-format position file
"""

import argparse


def read_gene_ids(fasta_path):
    """Gene ids = first token of each FASTA header."""
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
    """
    Aggregate each gene's span from its transcript lines: gene_id = mRNA Parent,
    span = (chromosome, min start, max end) over that gene's transcripts. Using
    transcript Parent keeps gene ids identical to the representative FASTA, which
    are derived the same way.
    """
    chromosome_of_gene = {}
    start_of_gene = {}
    end_of_gene = {}
    transcript_feature_types = {"mRNA", "transcript"}
    with open(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9 or columns[2] not in transcript_feature_types:
                continue
            attributes = parse_gff_attributes(columns[8])
            gene_id = attributes.get("Parent")
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
    parser.add_argument("--chrom-tag", required=True)
    parser.add_argument("--gene-prefix", default="")
    parser.add_argument("--out", required=True)
    arguments = parser.parse_args()

    genes_to_include = read_gene_ids(arguments.genes_from)
    chromosome_of_gene, start_of_gene, end_of_gene = read_gene_positions(arguments.gff)

    # give each distinct chromosome a stable integer index, in first-seen order
    chromosome_index = {}
    def tagged_chromosome(chromosome_name):
        if chromosome_name not in chromosome_index:
            chromosome_index[chromosome_name] = len(chromosome_index) + 1
        return f"{arguments.chrom_tag}{chromosome_index[chromosome_name]}"

    number_written = 0
    number_missing_position = 0
    with open(arguments.out, "w") as out_handle:
        for gene_id in genes_to_include:
            if gene_id not in chromosome_of_gene:
                number_missing_position += 1
                continue
            out_handle.write(
                f"{tagged_chromosome(chromosome_of_gene[gene_id])}\t"
                f"{arguments.gene_prefix}{gene_id}\t"
                f"{start_of_gene[gene_id]}\t{end_of_gene[gene_id]}\n"
            )
            number_written += 1

    assert number_written > 0, "no genes written to position file"
    print(f"  wrote {number_written} gene positions "
          f"({len(chromosome_index)} chromosomes, tag '{arguments.chrom_tag}')")
    if number_missing_position:
        print(f"  WARNING: {number_missing_position} proteins had no GFF position")


if __name__ == "__main__":
    main()
