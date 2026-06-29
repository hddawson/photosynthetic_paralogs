#!/usr/bin/env python3
"""
summarize_dupgen_types.py

DupGen_finder writes one file of duplicate PAIRS per mechanism
(<target>.wgd.pairs, .tandem.pairs, .proximal.pairs, .transposed.pairs,
.dispersed.pairs). A gene can appear in several. This collapses them to one
duplication type per gene using DupGen_finder's own priority order:

    wgd > tandem > proximal > transposed > dispersed

(i.e. a gene that is both a WGD and a dispersed duplicate is called WGD).

Output: gene_id <TAB> duplication_type, for every gene DupGen classified.
Genes not listed here were not classified as duplicated by DupGen (they are
singletons or unplaced) and are filled in as 'none' at the R join step.

Inputs:
  --dupgen-out-dir   directory containing <target>.<type>.pairs files
  --target           target species name (the file-name prefix)
Output:
  --out              per-gene duplication type table
"""

import argparse
import os

# highest priority first; this is the order DupGen_finder itself uses
DUPLICATION_TYPES_IN_PRIORITY_ORDER = [
    "wgd", "tandem", "proximal", "transposed", "dispersed"
]


def read_genes_in_pairs_file(pairs_file_path):
    """
    Return the set of target gene ids appearing in a .pairs file. The format is
    a header line then rows whose 1st and 3rd tab-columns are the two gene ids.
    Outgroup genes (prefixed OG_) are ignored -- only target genes are classified.
    """
    genes = set()
    if not os.path.exists(pairs_file_path):
        return genes
    with open(pairs_file_path) as handle:
        next(handle, None)  # skip header
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3:
                continue
            for gene_id in (columns[0], columns[2]):
                if gene_id and not gene_id.startswith("OG_"):
                    genes.add(gene_id)
    return genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dupgen-out-dir", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--out", required=True)
    arguments = parser.parse_args()

    # assign each gene the HIGHEST-priority type it appears in; because we iterate
    # from highest to lowest priority and never overwrite, first assignment wins
    duplication_type_of_gene = {}
    any_pairs_file_found = False
    for duplication_type in DUPLICATION_TYPES_IN_PRIORITY_ORDER:
        pairs_file_path = os.path.join(
            arguments.dupgen_out_dir, f"{arguments.target}.{duplication_type}.pairs"
        )
        if os.path.exists(pairs_file_path):
            any_pairs_file_found = True
        for gene_id in read_genes_in_pairs_file(pairs_file_path):
            duplication_type_of_gene.setdefault(gene_id, duplication_type)

    assert any_pairs_file_found, (
        f"no .pairs files found for target '{arguments.target}' in "
        f"{arguments.dupgen_out_dir} -- did DupGen_finder run?"
    )

    with open(arguments.out, "w") as out_handle:
        out_handle.write("gene_id\tduplication_type\n")
        for gene_id, duplication_type in duplication_type_of_gene.items():
            out_handle.write(f"{gene_id}\t{duplication_type}\n")

    # quick tally so the per-species counts are visible in the run log
    tally = {}
    for duplication_type in duplication_type_of_gene.values():
        tally[duplication_type] = tally.get(duplication_type, 0) + 1
    print(f"  classified {len(duplication_type_of_gene)} genes: " +
          ", ".join(f"{t}={tally.get(t, 0)}"
                    for t in DUPLICATION_TYPES_IN_PRIORITY_ORDER))


if __name__ == "__main__":
    main()
