#!/usr/bin/env python3
"""
Plot one hit-density ideogram per genome, with one strip per gene category.

Each contig row holds N stacked strips of binned density — one per category
present in that genome (the "null" housekeeping strip on top, then each
structure/PDB category). Each contig is divided into equal-width bins; each
bin cell is coloured by how many hits of that category fall inside it.
All strips share ONE colour scale and one colorbar, so brightness is directly
comparable across categories.

OUTPUT:
    results/gene_hits/genome_hit_density/<accession>_hit_density.png
"""

import zipfile
import os
import csv
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------- configuration ------------------------------- #

genome_zip_directory  = "data/hifi_angiosperm_genomes"
hits_summary_path     = "results/gene_hits/hits_summary.tsv"
gene_categories_path  = "data/reference_protein_sequences/gene_categories.tsv"
assembly_metadata_path = (
    "data/hifi_angiosperm_metadata/hifi_representative_newest_per_species.csv"
)
output_plot_directory = "results/gene_hits/genome_hit_density"

minimum_contig_length_nt = 1_000_000
maximum_contigs_to_plot  = 20
number_of_bins           = 10

# Single shared colormap for all strata so brightness is directly comparable.
density_colormap = "viridis"

# Fraction of each contig's row used by the stack of strips (rest is gap).
total_stack_height = 0.85

# --------------------------------------------------------------------------- #


def find_genome_fasta_inside_zip(zip_path):
    with zipfile.ZipFile(zip_path) as genome_zip:
        fna_entries = [
            entry for entry in genome_zip.namelist()
            if entry.endswith(".fna") and "genomic" in entry
        ]
        assert len(fna_entries) >= 1, f"No *genomic.fna in {zip_path}"
        if len(fna_entries) > 1:
            fna_entries = sorted(
                fna_entries,
                key=lambda entry: genome_zip.getinfo(entry).file_size,
                reverse=True,
            )
    return fna_entries[0]


def contig_lengths_from_zip(zip_path):
    fna_entry = find_genome_fasta_inside_zip(zip_path)
    contig_length_by_id = {}
    current_contig_id   = None
    current_length      = 0
    with zipfile.ZipFile(zip_path) as genome_zip:
        with genome_zip.open(fna_entry) as fasta_file:
            for raw_line in fasta_file:
                line = raw_line.decode("utf-8").rstrip()
                if line.startswith(">"):
                    if current_contig_id is not None:
                        contig_length_by_id[current_contig_id] = current_length
                    current_contig_id = line[1:].split()[0]
                    current_length    = 0
                else:
                    current_length += len(line)
    if current_contig_id is not None:
        contig_length_by_id[current_contig_id] = current_length
    return contig_length_by_id


def read_gene_categories(gene_categories_path):
    category_by_gene = {}
    with open(gene_categories_path) as tsv_file:
        for row in csv.DictReader(tsv_file, delimiter="\t"):
            category_by_gene[row["gene_name"]] = row["category"]
    return category_by_gene


def read_species_by_accession(assembly_metadata_path):
    species_by_accession = {}
    if not os.path.exists(assembly_metadata_path):
        return species_by_accession
    with open(assembly_metadata_path) as csv_file:
        for row in csv.DictReader(csv_file):
            accession = row.get("Assembly Accession", "").strip()
            organism  = row.get("Organism Name", "").strip()
            if accession and organism:
                species_by_accession[accession] = organism
                species_by_accession[accession.rsplit(".", 1)[0]] = organism
    return species_by_accession


def sanitise_gene_name(raw_gene_name):
    """
    Match the sanitisation applied by categorize_reference_genes.py:
    replace "|" and spaces with "_" so hits_summary gene names join to
    gene_categories gene names correctly.
    e.g. "sp|P00875|RBL_SPIOL" -> "sp_P00875_RBL_SPIOL"
    """
    return raw_gene_name.replace("|", "_").replace(" ", "_")


def read_hits_grouped_by_assembly(hits_summary_path, category_by_gene):
    """
    Return {assembly_accession: list of hit dicts}. Each hit carries its actual
    gene category in 'stratum' (genes with no category become 'uncategorised').
    Gene names are sanitised before the category lookup so pipe characters in
    hits_summary.tsv match underscore-sanitised names in gene_categories.tsv.
    """
    hits_by_assembly = {}
    with open(hits_summary_path) as tsv_file:
        for row in csv.DictReader(tsv_file, delimiter="\t"):
            hit_midpoint_nt = (int(row["genomic_start"]) + int(row["genomic_end"])) / 2
            sanitised_name  = sanitise_gene_name(row["gene_name"])
            category        = category_by_gene.get(sanitised_name, "uncategorised")
            hit = {
                "contig":          row["contig"],
                "hit_midpoint_nt": hit_midpoint_nt,
                "stratum":         category,
            }
            hits_by_assembly.setdefault(row["assembly_accession"], []).append(hit)
    return hits_by_assembly


def order_strata(strata_present):
    """
    Order strata for stacking: the housekeeping/null category first (top),
    then the rest sorted alphabetically.
    We treat both 'null' and any category containing 'housekeeping' as the
    top stratum, since the categoriser may label it from the filename.
    """
    null_strata  = sorted(s for s in strata_present
                          if s == "null" or "housekeeping" in s.lower())
    other_strata = sorted(s for s in strata_present
                          if s != "null" and "housekeeping" not in s.lower())
    return null_strata + other_strata


def count_hits_per_bin(contig_length, hit_positions, number_of_bins):
    bin_counts, _ = np.histogram(
        hit_positions, bins=number_of_bins, range=(0, contig_length)
    )
    return bin_counts


def plot_one_genome(assembly_accession, species_name, contig_length_by_id,
                    genome_hits, output_path):
    """Draw one density strip per category per contig, on a shared scale."""

    long_contigs = sorted(
        [(c, length) for c, length in contig_length_by_id.items()
         if length >= minimum_contig_length_nt],
        key=lambda pair: pair[1],
        reverse=True,
    )
    if not long_contigs:
        return False

    hidden_contig_count = max(0, len(long_contigs) - maximum_contigs_to_plot)
    plotted_contigs = long_contigs[:maximum_contigs_to_plot]

    # Which categories appear in this genome, ordered (null on top).
    strata_present = order_strata({hit["stratum"] for hit in genome_hits})
    n_strata = len(strata_present)
    if n_strata == 0:
        return False

    # Split hit positions by contig and category.
    positions_by_contig = {}
    for hit in genome_hits:
        positions_by_contig.setdefault(hit["contig"], {})
        positions_by_contig[hit["contig"]].setdefault(hit["stratum"], []).append(
            hit["hit_midpoint_nt"]
        )

    # Per-bin counts for every contig/stratum, plus one global max across all
    # strata so the single shared colour scale is meaningful.
    bin_counts = {}
    global_max_count = 1
    for contig_id, contig_length in plotted_contigs:
        for stratum in strata_present:
            positions = positions_by_contig.get(contig_id, {}).get(stratum, [])
            counts = count_hits_per_bin(contig_length, positions, number_of_bins)
            bin_counts[(contig_id, stratum)] = counts
            global_max_count = max(global_max_count, counts.max())

    # Even vertical spacing: N strips stacked within total_stack_height,
    # first stratum at the top of the row.
    strip_height = total_stack_height / n_strata
    strip_centre_offset = {
        stratum: total_stack_height / 2 - strip_height / 2 - index * strip_height
        for index, stratum in enumerate(strata_present)
    }

    colormap = plt.get_cmap(density_colormap)

    figure_height = max(3, len(plotted_contigs) * 0.25 * max(2, n_strata))
    figure, axis = plt.subplots(figsize=(11, figure_height))

    for row_index, (contig_id, contig_length) in enumerate(reversed(plotted_contigs)):
        bin_width_nt = contig_length / number_of_bins
        for stratum in strata_present:
            counts = bin_counts[(contig_id, stratum)]
            y_centre = row_index + strip_centre_offset[stratum]
            for bin_index in range(number_of_bins):
                colour_value = counts[bin_index] / global_max_count
                axis.barh(
                    y=y_centre,
                    width=bin_width_nt,
                    left=bin_index * bin_width_nt,
                    height=strip_height,
                    color=colormap(colour_value),
                    edgecolor="white",
                    linewidth=0.3,
                )
            # Short category tag at the left edge of each strip.
            axis.annotate(
                stratum[:8],
                xy=(0, y_centre), xytext=(-6, 0), textcoords="offset points",
                ha="right", va="center", fontsize=5, color="dimgrey",
            )

    axis.set_yticks(range(len(plotted_contigs)))
    axis.set_yticklabels([c for c, _ in reversed(plotted_contigs)])
    axis.set_xlabel("Position (nt)")
    axis.set_ylabel("Contig")
    axis.set_ylim(-1, len(plotted_contigs))

    title = assembly_accession
    if species_name:
        title = f"{species_name}  ({assembly_accession})"
    axis.set_title(
        f"{title} — hit density per {number_of_bins} bins\n"
        f"one strip per category (null on top): {', '.join(strata_present)}"
    )

    if hidden_contig_count > 0:
        axis.annotate(
            f"+ {hidden_contig_count} shorter contig(s) not shown",
            xy=(0, 0), xycoords="axes fraction",
            xytext=(0, -35), textcoords="offset points",
            ha="left", va="top", fontsize=8, style="italic", color="dimgrey",
        )

    scalar_mappable = plt.cm.ScalarMappable(
        cmap=colormap, norm=plt.Normalize(vmin=0, vmax=global_max_count)
    )
    colorbar = figure.colorbar(scalar_mappable, ax=axis, pad=0.01)
    colorbar.set_label("Hits per bin")

    figure.tight_layout()
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return True


def main():
    os.makedirs(output_plot_directory, exist_ok=True)

    category_by_gene     = read_gene_categories(gene_categories_path)
    species_by_accession = read_species_by_accession(assembly_metadata_path)
    hits_by_assembly     = read_hits_grouped_by_assembly(hits_summary_path, category_by_gene)
    assert len(hits_by_assembly) > 0, "No hits found in hits_summary.tsv"

    genome_zip_paths = sorted(glob.glob(os.path.join(genome_zip_directory, "*.zip")))
    assert len(genome_zip_paths) > 0, f"No .zip files in {genome_zip_directory}"
    zip_path_by_accession = {
        os.path.basename(zip_path).replace(".zip", ""): zip_path
        for zip_path in genome_zip_paths
    }

    plots_written = 0
    for assembly_accession, genome_hits in hits_by_assembly.items():
        zip_path = zip_path_by_accession.get(assembly_accession)
        if zip_path is None:
            print(f"  WARNING: no zip found for {assembly_accession}, skipping.")
            continue

        print(f"Reading contigs and plotting {assembly_accession}...")
        contig_length_by_id = contig_lengths_from_zip(zip_path)
        species_name = species_by_accession.get(assembly_accession, "")

        output_path = os.path.join(
            output_plot_directory, f"{assembly_accession}_hit_density.png"
        )
        if plot_one_genome(assembly_accession, species_name, contig_length_by_id,
                           genome_hits, output_path):
            plots_written += 1

    print(f"\nDone. Wrote {plots_written} plots to {output_plot_directory}")


if __name__ == "__main__":
    main()