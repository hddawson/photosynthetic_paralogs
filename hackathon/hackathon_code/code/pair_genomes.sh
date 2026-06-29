#!/usr/bin/env bash
#
# pair_genomes.sh
#
# For each GFF, find its genome FASTA and report whether it exists. Run ONCE
# before batching to see which species are runnable.
#
# Two GFF naming styles exist in this dataset and are handled here:
#   * snake_case full name : genus_species_GeneCAD_final.gff   (has underscores)
#         -> FASTA stem is <Genus-initial><species>            e.g. Athaliana
#   * already-abbreviated  : Carabica_GeneCAD_final.gff        (no underscore)
#         -> FASTA stem IS the GFF stem                        e.g. Carabica
#
# Genomes may be gzipped (.fa.gz) or plain (.fa) -- both are checked.
#
# Prints a tab-separated table to stdout:
#   gff_file <TAB> genome_path <TAB> status      (status = OK or MISSING)
# OK rows give a full genome path you can feed straight to find_duplicates.sh.
#
# Usage:
#   ./pair_genomes.sh <gff_dir> <genome_dir>

set -euo pipefail

gff_directory="${1:?need gff directory}"
genome_directory="${2:?need genome directory}"

printf "gff_file\tgenome_path\tstatus\n"
number_missing=0

for gff_path in "$gff_directory"/*_GeneCAD_final.gff; do
  gff_basename=$(basename "$gff_path")
  gff_stem=${gff_basename%_GeneCAD_final.gff}

  # decide the FASTA stem based on which GFF naming style this is
  if [[ "$gff_stem" == *_* ]]; then
    # snake_case: derive <Genus-initial><species>
    genus=${gff_stem%%_*}
    remainder=${gff_stem#*_}
    species=${remainder%%_*}
    genus_initial=${genus:0:1}; genus_initial=${genus_initial^^}
    fasta_stem="${genus_initial}${species}"
  else
    # already abbreviated: the GFF stem is the FASTA stem
    fasta_stem="$gff_stem"
  fi

  # look for the genome under either extension
  gzipped_candidate="${genome_directory}/${fasta_stem}.softmasked_all.fa.gz"
  plain_candidate="${genome_directory}/${fasta_stem}.softmasked_all.fa"
  if [ -e "$gzipped_candidate" ]; then
    printf "%s\t%s\tOK\n" "$gff_basename" "$gzipped_candidate"
  elif [ -e "$plain_candidate" ]; then
    printf "%s\t%s\tOK\n" "$gff_basename" "$plain_candidate"
  else
    printf "%s\t%s\tMISSING\n" "$gff_basename" "${fasta_stem}.softmasked_all.fa[.gz]"
    number_missing=$((number_missing + 1))
  fi
done

echo "" >&2
echo "genomes that did not resolve: ${number_missing}" >&2