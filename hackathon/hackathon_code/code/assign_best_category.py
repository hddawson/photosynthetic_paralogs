#!/usr/bin/env python3
"""
assign_best_category.py

Assign each genome gene a category (photosynthetic complex or housekeeping)
from its single best hit against the reference protein set.

The category of a reference protein comes from the user's gene_categories.tsv
(column: gene_name -> category). DIAMOND subject ids are sanitised the same way
that table was built (first token, '|' -> '_') so they join.

A genome gene is labelled only if its best hit clears the thresholds; genes with
no good hit are simply absent from the output (treated as 'unlabelled' later).

Inputs:
  --hits        DIAMOND outfmt 6: qseqid sseqid pident length evalue bitscore qcovhsp scovhsp
  --categories  gene_categories.tsv (columns include gene_name, category)
Output:
  --out         gene_id, category, reference_hit, evalue, percent_identity, query_coverage
"""

import argparse

# a genome gene is labelled only if its best hit passes all of these
MAXIMUM_EVALUE = 1e-10
MINIMUM_QUERY_COVERAGE_PERCENT = 50
MINIMUM_PERCENT_IDENTITY = 30


def sanitise_reference_id(raw_subject_id):
    """Match the reference-table convention: first token, '|' and ' ' -> '_'."""
    first_token = raw_subject_id.split()[0]
    return first_token.replace("|", "_").replace(" ", "_")


def read_category_of_reference(categories_path):
    """Build {sanitised_reference_gene_name: category} from gene_categories.tsv."""
    category_of_reference = {}
    with open(categories_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        name_column = header.index("gene_name")
        category_column = header.index("category")
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) <= max(name_column, category_column):
                continue
            category_of_reference[columns[name_column]] = columns[category_column]
    assert category_of_reference, f"no categories read from {categories_path}"
    return category_of_reference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True)
    parser.add_argument("--categories", required=True)
    parser.add_argument("--out", required=True)
    arguments = parser.parse_args()

    category_of_reference = read_category_of_reference(arguments.categories)

    # for each genome gene keep its best passing hit (lowest e-value, then highest
    # bitscore); store (evalue, bitscore, reference_id, pident, qcov)
    best_hit_for_gene = {}
    with open(arguments.hits) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            query_gene = fields[0]
            reference_id = fields[1]
            percent_identity = float(fields[2])
            evalue = float(fields[4])
            bitscore = float(fields[5])
            query_coverage = float(fields[6])

            if evalue > MAXIMUM_EVALUE:
                continue
            if query_coverage < MINIMUM_QUERY_COVERAGE_PERCENT:
                continue
            if percent_identity < MINIMUM_PERCENT_IDENTITY:
                continue

            candidate = (evalue, -bitscore, reference_id, percent_identity, query_coverage)
            # smaller evalue wins; for ties, larger bitscore (stored negated) wins
            if query_gene not in best_hit_for_gene or candidate < best_hit_for_gene[query_gene]:
                best_hit_for_gene[query_gene] = candidate

    number_with_unknown_reference = 0
    with open(arguments.out, "w") as out_handle:
        out_handle.write("gene_id\tcategory\treference_hit\tevalue\t"
                         "percent_identity\tquery_coverage\n")
        for gene_id, (evalue, negative_bitscore, reference_id,
                      percent_identity, query_coverage) in best_hit_for_gene.items():
            reference_name = sanitise_reference_id(reference_id)
            category = category_of_reference.get(reference_name)
            if category is None:
                number_with_unknown_reference += 1
                category = "unknown_reference"
            out_handle.write(
                f"{gene_id}\t{category}\t{reference_name}\t{evalue:.3g}\t"
                f"{percent_identity:.1f}\t{query_coverage:.1f}\n"
            )

    print(f"  labelled {len(best_hit_for_gene)} genes")
    if number_with_unknown_reference:
        print(f"  WARNING: {number_with_unknown_reference} best hits had no category "
              f"in the reference table")


if __name__ == "__main__":
    main()
