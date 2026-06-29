#!/usr/bin/env python3
"""
cluster_into_families.py

Group genes into families from an all-vs-all DIAMOND search, then emit a
per-gene duplication table.

Method: build a graph where each gene is a node and an edge is drawn between two
DIFFERENT genes whose alignment passes the thresholds below. Each connected
component of that graph is one gene family. A gene is "duplicated" if its family
has more than one member.

NOTE on method: connected components = single-linkage clustering. It is simple
and dependency-free, but it can "chain" (A-B, B-C links A and C even if A and C
are unrelated). The coverage + identity thresholds below are what keep chaining
in check; loosen them and families merge, tighten them and families split. If
you later want stricter families, swap this step for MCL or MMseqs2 clustering.

Inputs:
  --hits        DIAMOND outfmt 6 with columns:
                qseqid sseqid pident length evalue bitscore qcovhsp scovhsp
  --gene-fasta  representative proteins (gives the FULL gene set, so singletons
                are included as families of size 1)
Outputs:
  --out-per-gene  gene_id, family_id, family_size, is_duplicated
  --out-families  family_id, family_size, member_gene_ids(comma-separated)

Standard library only.
"""

import argparse

# ---- linking thresholds (edit here) -----------------------------------------
# An edge is drawn between two genes only if a hit clears ALL of these.
MAXIMUM_EVALUE = 1e-5          # significance of the alignment
MINIMUM_PERCENT_IDENTITY = 30  # %, a common homology floor for paralog families
MINIMUM_COVERAGE_PERCENT = 50  # % of BOTH query and subject must align,
                               # so a single shared domain isn't enough to link


def read_gene_ids_from_fasta(fasta_path):
    """Collect every gene id (FASTA header first token) so singletons are kept."""
    gene_ids = []
    with open(fasta_path) as fasta_handle:
        for line in fasta_handle:
            if line.startswith(">"):
                gene_ids.append(line[1:].split()[0])
    assert len(gene_ids) > 0, "no gene ids in representative FASTA"
    return gene_ids


class UnionFind:
    """Minimal union-find (disjoint set) for connected components."""

    def __init__(self, items):
        # each item starts as its own parent (its own component)
        self.parent_of = {item: item for item in items}

    def find_root(self, item):
        # walk to the component root, compressing the path as we go
        root = item
        while self.parent_of[root] != root:
            root = self.parent_of[root]
        while self.parent_of[item] != root:
            self.parent_of[item], item = root, self.parent_of[item]
        return root

    def union(self, first_item, second_item):
        self.parent_of[self.find_root(first_item)] = self.find_root(second_item)


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--hits", required=True)
    argument_parser.add_argument("--gene-fasta", required=True)
    argument_parser.add_argument("--out-per-gene", required=True)
    argument_parser.add_argument("--out-families", required=True)
    arguments = argument_parser.parse_args()

    all_gene_ids = read_gene_ids_from_fasta(arguments.gene_fasta)
    components = UnionFind(all_gene_ids)

    # draw an edge for every hit that clears all thresholds and is not a self-hit
    number_of_linking_hits = 0
    with open(arguments.hits) as hits_handle:
        for line in hits_handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            query_gene, subject_gene = fields[0], fields[1]
            percent_identity = float(fields[2])
            evalue = float(fields[4])
            query_coverage = float(fields[6])
            subject_coverage = float(fields[7])

            if query_gene == subject_gene:
                continue
            if evalue > MAXIMUM_EVALUE:
                continue
            if percent_identity < MINIMUM_PERCENT_IDENTITY:
                continue
            if query_coverage < MINIMUM_COVERAGE_PERCENT:
                continue
            if subject_coverage < MINIMUM_COVERAGE_PERCENT:
                continue

            components.union(query_gene, subject_gene)
            number_of_linking_hits += 1

    # gather genes by their component root = family
    genes_in_family = {}
    for gene_id in all_gene_ids:
        family_root = components.find_root(gene_id)
        genes_in_family.setdefault(family_root, []).append(gene_id)

    # give families stable numeric ids, largest first
    families_sorted_by_size = sorted(
        genes_in_family.values(), key=len, reverse=True
    )

    # sanity check: every gene lands in exactly one family
    total_genes_assigned = sum(len(members) for members in families_sorted_by_size)
    assert total_genes_assigned == len(all_gene_ids), (
        f"gene count mismatch: {total_genes_assigned} assigned vs "
        f"{len(all_gene_ids)} input genes"
    )

    family_id_of_gene = {}
    family_size_of_gene = {}
    with open(arguments.out_families, "w") as families_output:
        families_output.write("family_id\tfamily_size\tmember_gene_ids\n")
        for family_index, member_gene_ids in enumerate(families_sorted_by_size, start=1):
            family_identifier = f"family_{family_index:06d}"
            family_size = len(member_gene_ids)
            for gene_id in member_gene_ids:
                family_id_of_gene[gene_id] = family_identifier
                family_size_of_gene[gene_id] = family_size
            families_output.write(
                f"{family_identifier}\t{family_size}\t{','.join(member_gene_ids)}\n"
            )

    # per-gene table: this is the file you join your category labels onto later
    number_of_duplicated_genes = 0
    with open(arguments.out_per_gene, "w") as per_gene_output:
        per_gene_output.write("gene_id\tfamily_id\tfamily_size\tis_duplicated\n")
        for gene_id in all_gene_ids:
            family_size = family_size_of_gene[gene_id]
            is_duplicated = family_size > 1
            number_of_duplicated_genes += int(is_duplicated)
            per_gene_output.write(
                f"{gene_id}\t{family_id_of_gene[gene_id]}\t{family_size}\t"
                f"{'TRUE' if is_duplicated else 'FALSE'}\n"
            )

    print(f"  linking hits kept            : {number_of_linking_hits}")
    print(f"  genes total                  : {len(all_gene_ids)}")
    print(f"  genes in a family of size >1 : {number_of_duplicated_genes}")
    print(f"  families total               : {len(families_sorted_by_size)}")


if __name__ == "__main__":
    main()
