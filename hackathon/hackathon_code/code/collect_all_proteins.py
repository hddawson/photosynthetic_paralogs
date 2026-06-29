#!/usr/bin/env python3
"""
collect_all_proteins.py

Concatenate the representative_proteins.faa from every successfully-run species
into one combined FASTA for transit-peptide prediction.

Each sequence header is prefixed with the species name so hits can be traced back:
  >Pvaginatum__Chr01_gene_4
  >arabidopsis_thaliana__AT1G01010.1
  ...

Usage:
    python3 code/collect_all_proteins.py [results_dir] [output.faa]
"""

import os
import sys

results_directory   = sys.argv[1] if len(sys.argv) > 1 else "results"
output_fasta_path   = sys.argv[2] if len(sys.argv) > 2 else "data/all_genecad_proteins.faa"

number_of_species  = 0
number_of_proteins = 0

with open(output_fasta_path, "w") as output_handle:
    for species_folder_name in sorted(os.listdir(results_directory)):
        species_folder_path    = os.path.join(results_directory, species_folder_name)
        representative_proteins_path = os.path.join(
            species_folder_path, "representative_proteins.faa"
        )
        if not os.path.isfile(representative_proteins_path):
            continue

        number_of_species += 1
        with open(representative_proteins_path) as input_handle:
            for line in input_handle:
                if line.startswith(">"):
                    # prefix the gene id with the species folder name
                    gene_id = line[1:].split()[0]
                    output_handle.write(f">{species_folder_name}__{gene_id}\n")
                    number_of_proteins += 1
                else:
                    output_handle.write(line)

assert number_of_proteins > 0, "no proteins written -- check results directory"
print(f"species collected : {number_of_species}")
print(f"proteins written  : {number_of_proteins}")
print(f"output            : {output_fasta_path}")
