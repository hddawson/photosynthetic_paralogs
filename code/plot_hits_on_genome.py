#!/usr/bin/env python3
"""
Plot one ideogram per genome: chromosomes as horizontal bars with miniprot
gene hits marked along them.

Hits are coloured by category (from gene_categories.tsv): the "null"
housekeeping genes share one colour; each structure category gets its own.
Contigs are capped to the longest N per genome to keep plots readable, with
a note of how many were hidden. Each plot is titled with the species name
pulled from the assembly metadata CSV.

OUTPUT:
    results/gene_hits/genome_hit_maps/<accession>_hit_map.png
"""

import zipfile
import os
import csv
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------- configuration ------------------------------- #

genome_zip_directory  = "data/hifi_angiosperm_genomes"
hits_summary_path     = "results/gene_hits/hits_summary.tsv"
gene_categories_path  = "data/reference_protein_sequences/gene_categories.tsv"
# Assembly metadata CSV with accession + organism name (from the NCBI step).
assembly_metadata_path = (
    "data/hifi_angiosperm_metadata/hifi_representative_newest_per_species.csv"
)
output_plot_directory = "results/gene_hits/genome_hit_maps"

# Drop contigs shorter than this (tiny unplaced scaffolds). 0 shows everything.
minimum_contig_length_nt = 1_000_000

# Show at most this many contigs (the longest ones) per genome.
maximum_contigs_to_plot = 20

# --------------------------------------------------------------------------- #


def find_genome_fasta_inside_zip(zip_path):
    """Return the largest *genomic.fna entry inside an NCBI datasets zip."""
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
    """Stream the genome FASTA from the zip; return {contig_id: length_nt}."""
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
    """Return {gene_name: category} from gene_categories.tsv."""
    category_by_gene = {}
    with open(gene_categories_path) as tsv_file:
        for row in csv.DictReader(tsv_file, delimiter="\t"):
            category_by_gene[row["gene_name"]] = row["category"]
    return category_by_gene


def read_species_by_accession(assembly_metadata_path):
    species_by_accession = {}
    if not os.path.exists(assembly_metadata_path):
        print(f"  NOTE: metadata CSV not found ({assembly_metadata_path})")
        return species_by_accession

    with open(assembly_metadata_path) as csv_file:
        for row in csv.DictReader(csv_file):
            accession = row.get("Assembly Accession", "").strip()
            organism  = row.get("Organism Name", "").strip()
            if accession and organism:
                # Store both the full versioned accession and the bare accession
                # so we match regardless of whether the hits carry a version suffix.
                species_by_accession[accession] = organism
                bare_accession = accession.rsplit(".", 1)[0]
                species_by_accession[bare_accession] = organism

    if not species_by_accession:
        print("  WARNING: no accession/organism pairs read from metadata CSV. "
              "Check column names are 'Assembly Accession' and 'Organism Name'.")
    return species_by_accession


def read_hits_grouped_by_assembly(hits_summary_path):
    """Return {assembly_accession: list of hit dicts}."""
    hits_by_assembly = {}
    with open(hits_summary_path) as tsv_file:
        for row in csv.DictReader(tsv_file, delimiter="\t"):
            hit_midpoint_nt = (int(row["genomic_start"]) + int(row["genomic_end"])) / 2
            hit = {
                "contig":          row["contig"],
                "hit_midpoint_nt": hit_midpoint_nt,
                "gene_name":       row["gene_name"],
            }
            hits_by_assembly.setdefault(row["assembly_accession"], []).append(hit)
    return hits_by_assembly


def build_category_colour_map(category_by_gene):
    """
    Assign a colour to each category. "null" is fixed grey; every other
    category (structure) gets a distinct colour from tab10.
    Returns {category: colour}.
    """
    non_null_categories = sorted(
        {c for c in category_by_gene.values() if c != "null"}
    )
    palette = plt.get_cmap("tab10")
    colour_by_category = {"null": "grey"}
    for index, category in enumerate(non_null_categories):
        colour_by_category[category] = palette(index % 10)
    return colour_by_category


def plot_one_genome(assembly_accession, species_name, contig_length_by_id,
                    genome_hits, category_by_gene, colour_by_category,
                    output_path):
    """Draw the ideogram for one genome and save a PNG. Returns True if drawn."""

    # Keep long-enough contigs, sorted longest-first.
    long_contigs = sorted(
        [(c, length) for c, length in contig_length_by_id.items()
         if length >= minimum_contig_length_nt],
        key=lambda pair: pair[1],
        reverse=True,
    )
    if not long_contigs:
        return False

    # Cap to the longest N contigs; remember how many we hid.
    hidden_contig_count = max(0, len(long_contigs) - maximum_contigs_to_plot)
    plotted_contigs = long_contigs[:maximum_contigs_to_plot]

    # Longest contig at the top of the plot.
    contig_y_position = {
        contig_id: row_index
        for row_index, (contig_id, _) in enumerate(reversed(plotted_contigs))
    }

    figure_height = max(3, len(plotted_contigs) * 0.3)
    figure, axis = plt.subplots(figsize=(10, figure_height))

    # Chromosome bars.
    for contig_id, contig_length in plotted_contigs:
        y = contig_y_position[contig_id]
        axis.plot([0, contig_length], [y, y],
                  color="lightgrey", linewidth=4, solid_capstyle="round", zorder=1)

    # Hit ticks, coloured by category. Genes with no category fall back to black.
    categories_present = set()
    for hit in genome_hits:
        if hit["contig"] not in contig_y_position:
            continue
        category = category_by_gene.get(hit["gene_name"], "uncategorised")
        categories_present.add(category)
        colour = colour_by_category.get(category, "black")
        y = contig_y_position[hit["contig"]]
        axis.scatter(hit["hit_midpoint_nt"], y,
                     marker="|", s=120, color=colour, zorder=2)

    axis.set_yticks(range(len(plotted_contigs)))
    axis.set_yticklabels([c for c, _ in reversed(plotted_contigs)])
    axis.set_xlabel("Position (nt)")
    axis.set_ylabel("Contig")

    # Title: species name (if known) plus the accession.
    title = assembly_accession
    if species_name:
        title = f"{species_name}  ({assembly_accession})"
    axis.set_title(title)

    # Note hidden contigs below the plot.
    if hidden_contig_count > 0:
        axis.annotate(
            f"+ {hidden_contig_count} shorter contig(s) not shown",
            xy=(0, 0), xycoords="axes fraction",
            xytext=(0, -35), textcoords="offset points",
            ha="left", va="top", fontsize=8, style="italic", color="dimgrey",
        )

    # Legend: one entry per category actually shown.
    legend_handles = [
        plt.Line2D([0], [0], marker="|", linestyle="None",
                   color=colour_by_category.get(category, "black"),
                   label=category, markersize=10)
        for category in sorted(categories_present)
    ]
    if legend_handles:
        axis.legend(handles=legend_handles, title="Category",
                    bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return True


def main():
    os.makedirs(output_plot_directory, exist_ok=True)

    category_by_gene     = read_gene_categories(gene_categories_path)
    species_by_accession = read_species_by_accession(assembly_metadata_path)
    hits_by_assembly     = read_hits_grouped_by_assembly(hits_summary_path)
    assert len(hits_by_assembly) > 0, "No hits found in hits_summary.tsv"

    colour_by_category = build_category_colour_map(category_by_gene)

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
            output_plot_directory, f"{assembly_accession}_hit_map.png"
        )
        if plot_one_genome(assembly_accession, species_name, contig_length_by_id,
                           genome_hits, category_by_gene, colour_by_category,
                           output_path):
            plots_written += 1

    print(f"\nDone. Wrote {plots_written} plots to {output_plot_directory}")


if __name__ == "__main__":
    main()