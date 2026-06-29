#!/usr/bin/env python3
"""
Search for copies of one or more genes across HiFi angiosperm genomes using
miniprot (protein-to-genome spliced alignment).

All protein FASTA files found in reference_proteins_directory are concatenated
into a single query file and passed to miniprot in one batch per genome.
This is fast because miniprot indexes the genome once regardless of how many
query proteins are used.

Hits are labelled by gene name (taken from the FASTA header of the reference
protein) and written to a single all_hits.fna output file. The summary TSV
has one row per hit, with a gene_name column so you can split by gene in R.

To add a new gene to the search, just drop its protein FASTA into
reference_proteins_directory and re-run.

REQUIREMENTS:
    miniprot at /programs/miniprot-0.13/miniprot  (BioHPC)
    Python standard library only.

OUTPUT (all in results_directory/):
    all_hits.fna       — genomic locus sequences for all genes, all genomes
    hits_summary.tsv   — one row per hit: assembly, gene, coordinates, etc.
"""

import subprocess
import zipfile
import os
import csv
import glob

# ----------------------------- configuration ------------------------------- #

# Genome zip files downloaded by download_hifi_angiosperm_genomes.py.
genome_zip_directory = "data/hifi_angiosperm_genomes"

# Directory containing one or more reference protein FASTA files (.fa / .faa / .fasta).
# All files in this directory are concatenated into a single query for miniprot.
reference_proteins_directory = "data/reference_protein_sequences"

# Where to write results.
results_directory = "results/gene_hits"

# Full path to miniprot binary (BioHPC).
miniprot_binary = "/programs/miniprot-0.13/miniprot"

# Number of CPU threads for miniprot.
miniprot_threads = 10

# Minimum fraction of the query protein that must align to retain a hit.
# E.g. 0.5 means at least half the reference protein must be covered.
minimum_query_coverage_fraction = 0.5

# Nucleotides of genomic context to extract on each side of each hit locus.
# This flanking sequence is what Tiberius uses to re-predict the gene model.
flanking_nucleotides = 1500

# --------------------------------------------------------------------------- #


def collect_and_concatenate_reference_proteins(reference_proteins_directory,
                                               output_fasta_path):
    """
    Find all FASTA files in reference_proteins_directory, concatenate them
    into a single file at output_fasta_path, and return a dict mapping
    each protein's sequence ID to its length in amino acids.

    The per-protein length dict is used later to compute query coverage
    fractions, which differ per protein.
    """
    fasta_extensions = ("*.fa", "*.faa", "*.fasta")
    reference_fasta_paths = []
    for extension in fasta_extensions:
        reference_fasta_paths.extend(
            glob.glob(os.path.join(reference_proteins_directory, extension))
        )
    reference_fasta_paths = sorted(reference_fasta_paths)

    assert len(reference_fasta_paths) > 0, (
        f"No FASTA files found in {reference_proteins_directory}. "
        "Add at least one .fa/.faa/.fasta reference protein file."
    )
    print(f"Found {len(reference_fasta_paths)} reference protein file(s):")
    for path in reference_fasta_paths:
        print(f"  {path}")

    # Concatenate all reference files and build the protein length lookup.
    protein_lengths_by_id = {}
    with open(output_fasta_path, "w") as combined_fasta:
        for fasta_path in reference_fasta_paths:
            current_protein_id = None
            current_length     = 0
            with open(fasta_path) as individual_fasta:
                for line in individual_fasta:
                    line = line.rstrip()
                    if line.startswith(">"):
                        # Save the previous protein's length before moving on.
                        if current_protein_id is not None:
                            protein_lengths_by_id[current_protein_id] = current_length
                        # Protein ID is the first whitespace-delimited token after ">".
                        current_protein_id = line[1:].split()[0]
                        current_length     = 0
                    else:
                        current_length += len(line)
                    combined_fasta.write(line + "\n")
                # Save the last protein in this file.
                if current_protein_id is not None:
                    protein_lengths_by_id[current_protein_id] = current_length

    assert len(protein_lengths_by_id) > 0, "No protein sequences were read."
    print(f"Total reference proteins: {len(protein_lengths_by_id)}")
    return protein_lengths_by_id


def find_genome_fasta_inside_zip(zip_path):
    """
    Return the zip-internal path of the largest *genomic.fna file inside an
    NCBI datasets genome zip. Largest by compressed size = primary assembly.
    """
    with zipfile.ZipFile(zip_path) as genome_zip:
        fna_entries = [
            entry for entry in genome_zip.namelist()
            if entry.endswith(".fna") and "genomic" in entry
        ]
    assert len(fna_entries) >= 1, (
        f"No *genomic.fna found in {zip_path}."
    )
    if len(fna_entries) > 1:
        with zipfile.ZipFile(zip_path) as genome_zip:
            fna_entries = sorted(
                fna_entries,
                key=lambda entry: genome_zip.getinfo(entry).file_size,
                reverse=True,
            )
    return fna_entries[0]


def extract_genome_fasta(zip_path, fna_entry_name, extract_to_directory):
    """Extract one .fna entry from the zip. Returns the path of the extracted file."""
    with zipfile.ZipFile(zip_path) as genome_zip:
        genome_zip.extract(fna_entry_name, extract_to_directory)
    extracted_path = os.path.join(extract_to_directory, fna_entry_name)
    assert os.path.exists(extracted_path), f"Extraction failed: {extracted_path}"
    return extracted_path


def run_miniprot(combined_query_fasta_path, genome_fasta_path,
                 gff_output_path, threads):
    """
    Run miniprot with all reference proteins against one genome in a single pass.

    --gff   : produce GFF3 output (one mRNA + CDS features per hit locus)
    --outs  : retain hits scoring >= 90% of the top hit per locus, so that
              paralog copies are not suppressed by the best-hit filter
    -t      : threads
    """
    with open(gff_output_path, "w") as gff_file:
        completed = subprocess.run(
            [
                miniprot_binary,
                "--gff",
                "--outs=0.9",
                "-t", str(threads),
                genome_fasta_path,
                combined_query_fasta_path,
            ],
            stdout=gff_file,
            stderr=subprocess.PIPE,
            text=True,
        )
    assert completed.returncode == 0, f"miniprot failed:\n{completed.stderr}"


def parse_miniprot_gff3(gff_path, protein_lengths_by_id):
    """
    Parse miniprot GFF3 output into a list of hit locus dicts.

    Each mRNA line in the GFF has a Target= attribute that records which
    reference protein matched and what portion of it was aligned:
        Target=<protein_id> <aligned_start_aa> <aligned_end_aa>

    We use this to identify the gene name and compute query coverage
    (aligned aa / total protein length) for each hit.

    Hits below minimum_query_coverage_fraction are discarded.
    CDS child features are collected per locus for exon counting.
    """
    loci_by_id = {}

    with open(gff_path) as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue

            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue

            contig, source, feature_type, start_str, end_str, \
                score_str, strand, phase, attributes = columns

            # Parse "KEY=VALUE;KEY=VALUE;..." into a dict.
            attribute_dict = {}
            for pair in attributes.split(";"):
                pair = pair.strip()
                if "=" in pair:
                    key, value = pair.split("=", 1)
                    attribute_dict[key] = value

            if feature_type == "mRNA":
                locus_id = attribute_dict.get("ID", f"locus_{len(loci_by_id)}")

                # Extract the matched reference protein ID and aligned range.
                # Target attribute format: "<protein_id> <aa_start> <aa_end>"
                gene_name               = "unknown"
                query_coverage_fraction = 0.0
                target_field = attribute_dict.get("Target", "")
                if target_field:
                    target_parts = target_field.split()
                    if len(target_parts) == 3:
                        gene_name    = target_parts[0]
                        aligned_aa   = int(target_parts[2]) - int(target_parts[1])
                        protein_length = protein_lengths_by_id.get(gene_name, 1)
                        query_coverage_fraction = aligned_aa / protein_length

                loci_by_id[locus_id] = {
                    "locus_id":                locus_id,
                    "gene_name":               gene_name,
                    "contig":                  contig,
                    "genomic_start":           int(start_str),
                    "genomic_end":             int(end_str),
                    "strand":                  strand,
                    "score":                   score_str,
                    "query_coverage_fraction": round(query_coverage_fraction, 3),
                    "cds_intervals":           [],
                }

            elif feature_type == "CDS":
                parent_id = attribute_dict.get("Parent")
                if parent_id and parent_id in loci_by_id:
                    loci_by_id[parent_id]["cds_intervals"].append(
                        (int(start_str), int(end_str))
                    )

    # Keep only loci meeting the coverage threshold; sort exons by position.
    passing_loci = []
    for locus in loci_by_id.values():
        if locus["query_coverage_fraction"] >= minimum_query_coverage_fraction:
            locus["cds_intervals"].sort(key=lambda interval: interval[0])
            passing_loci.append(locus)

    return passing_loci


def load_genome_sequences(genome_fasta_path):
    """
    Read a genome FASTA file into a dict: {contig_id: nucleotide_sequence}.
    Contig ID is the first whitespace-delimited token of each header line.
    """
    sequences     = {}
    current_id    = None
    current_parts = []

    with open(genome_fasta_path) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_parts)
                current_id    = line[1:].split()[0]
                current_parts = []
            elif line:
                current_parts.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_parts)

    assert len(sequences) > 0, f"No sequences loaded from {genome_fasta_path}"
    return sequences


def extract_genomic_locus_with_flanking(locus, genome_sequences,
                                         assembly_accession, copy_number):
    """
    Extract the full genomic span of one hit locus (including introns) plus
    flanking_nucleotides of context on each side.

    GFF3 uses 1-based inclusive coordinates; Python slicing is 0-based exclusive,
    so we subtract 1 from the start and clamp both ends to the contig length.

    The flanking context is included so Tiberius can see splice sites and
    upstream/downstream regulatory sequence when re-predicting the gene model.

    Header format: {accession}_{gene_name}_copy{N} contig=... locus=... etc.
    Returns (fasta_header, sequence) or None if the contig is not in the genome.
    """
    contig_sequence = genome_sequences.get(locus["contig"])
    if contig_sequence is None:
        print(f"    WARNING: contig {locus['contig']} not found — skipping.")
        return None

    contig_length = len(contig_sequence)

    # Convert from 1-based GFF3 to 0-based Python and add flanking context.
    extract_start = max(0, locus["genomic_start"] - 1 - flanking_nucleotides)
    extract_end   = min(contig_length, locus["genomic_end"] + flanking_nucleotides)

    extracted_sequence = contig_sequence[extract_start:extract_end]

    # Sanitise the gene name for use in the FASTA header (| and spaces cause
    # problems in downstream tools like RAxML).
    safe_gene_name = locus["gene_name"].replace("|", "_").replace(" ", "_")

    fasta_header = (
        f"{assembly_accession}_{safe_gene_name}_copy{copy_number} "
        f"gene={locus['gene_name']} "
        f"contig={locus['contig']} "
        f"locus={locus['genomic_start']}-{locus['genomic_end']} "
        f"strand={locus['strand']} "
        f"extracted={extract_start + 1}-{extract_end} "
        f"query_coverage={locus['query_coverage_fraction']:.2f} "
        f"n_exons={len(locus['cds_intervals'])}"
    )

    return fasta_header, extracted_sequence


def process_one_genome(zip_path, combined_query_fasta_path, protein_lengths_by_id,
                       working_directory, all_hits_fasta_handle):
    """
    Run the full pipeline for one genome zip file:
        extract → miniprot → parse GFF → extract sequences → write to shared FASTA.

    Returns a list of summary row dicts (one per passing hit locus).
    Sequences are appended directly to all_hits_fasta_handle so we only keep
    one genome in memory at a time.
    """
    assembly_accession = os.path.basename(zip_path).replace(".zip", "")
    genome_work_dir    = os.path.join(working_directory, assembly_accession)
    os.makedirs(genome_work_dir, exist_ok=True)

    print(f"\n  [{assembly_accession}] Extracting genome FASTA...")
    fna_entry         = find_genome_fasta_inside_zip(zip_path)
    genome_fasta_path = extract_genome_fasta(zip_path, fna_entry, genome_work_dir)

    print(f"  [{assembly_accession}] Running miniprot...")
    gff_output_path = os.path.join(genome_work_dir, "miniprot_hits.gff")
    run_miniprot(combined_query_fasta_path, genome_fasta_path,
                 gff_output_path, miniprot_threads)

    passing_loci = parse_miniprot_gff3(gff_output_path, protein_lengths_by_id)
    print(f"  [{assembly_accession}] Loci passing coverage filter: {len(passing_loci)}")

    if not passing_loci:
        return []

    genome_sequences = load_genome_sequences(genome_fasta_path)

    # Count copies per gene independently so numbering resets per gene.
    copy_counter_by_gene = {}
    summary_rows         = []

    for locus in passing_loci:
        gene_name   = locus["gene_name"]
        copy_number = copy_counter_by_gene.get(gene_name, 0) + 1
        copy_counter_by_gene[gene_name] = copy_number

        result = extract_genomic_locus_with_flanking(
            locus, genome_sequences, assembly_accession, copy_number
        )
        if result is None:
            continue

        fasta_header, sequence = result
        all_hits_fasta_handle.write(f">{fasta_header}\n{sequence}\n")

        summary_rows.append({
            "assembly_accession":      assembly_accession,
            "gene_name":               gene_name,
            "copy_number":             copy_number,
            "contig":                  locus["contig"],
            "genomic_start":           locus["genomic_start"],
            "genomic_end":             locus["genomic_end"],
            "strand":                  locus["strand"],
            "locus_length_nt":         locus["genomic_end"] - locus["genomic_start"] + 1,
            "n_exons":                 len(locus["cds_intervals"]),
            "query_coverage_fraction": locus["query_coverage_fraction"],
            "miniprot_score":          locus["score"],
        })

    return summary_rows


def main():
    assert os.path.exists(miniprot_binary), \
        f"miniprot not found at {miniprot_binary}."
    assert os.path.isdir(reference_proteins_directory), \
        f"Reference proteins directory not found: {reference_proteins_directory}"

    os.makedirs(results_directory, exist_ok=True)
    working_directory = os.path.join(results_directory, "working")
    os.makedirs(working_directory, exist_ok=True)

    # Concatenate all reference protein FASTAs into one combined query file.
    combined_query_fasta_path = os.path.join(working_directory, "all_reference_proteins.faa")
    protein_lengths_by_id = collect_and_concatenate_reference_proteins(
        reference_proteins_directory, combined_query_fasta_path
    )

    genome_zip_paths = sorted(glob.glob(os.path.join(genome_zip_directory, "*.zip")))
    assert len(genome_zip_paths) > 0, \
        f"No .zip files found in {genome_zip_directory}."
    print(f"\nFound {len(genome_zip_paths)} genome zip files to process.")

    all_hits_fasta_path = os.path.join(results_directory, "all_hits.fna")
    summary_tsv_path    = os.path.join(results_directory, "hits_summary.tsv")

    all_summary_rows = []

    # Open the shared output FASTA once and pass the handle into each genome's
    # processing function so all hits are appended in genome order.
    with open(all_hits_fasta_path, "w") as all_hits_fasta_handle:
        for zip_path in genome_zip_paths:
            print(f"\nProcessing: {zip_path}")
            summary_rows = process_one_genome(
                zip_path, combined_query_fasta_path, protein_lengths_by_id,
                working_directory, all_hits_fasta_handle,
            )
            all_summary_rows.extend(summary_rows)

    tsv_column_names = [
        "assembly_accession", "gene_name", "copy_number", "contig",
        "genomic_start", "genomic_end", "strand",
        "locus_length_nt", "n_exons",
        "query_coverage_fraction", "miniprot_score",
    ]
    with open(summary_tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_column_names, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_summary_rows)

    total_hits        = len(all_summary_rows)
    genomes_with_hits = len({r["assembly_accession"] for r in all_summary_rows})
    genes_found       = sorted({r["gene_name"] for r in all_summary_rows})

    print(f"\nDone.")
    print(f"  Genomes with at least one hit : {genomes_with_hits} / {len(genome_zip_paths)}")
    print(f"  Total hit loci                : {total_hits}")
    print(f"  Genes recovered               : {genes_found}")
    print(f"  All hit sequences             : {all_hits_fasta_path}")
    print(f"  Summary TSV                   : {summary_tsv_path}")


if __name__ == "__main__":
    main()