#!/usr/bin/env python3
"""
Find copies of rbcS (Rubisco small subunit) in HiFi angiosperm genomes
using miniprot (protein-to-genome spliced alignment).

WHY miniprot OVER tBLASTn:
    miniprot is purpose-built for aligning a protein query to a genomic
    sequence. It handles introns natively, so each hit it reports is already
    a complete spliced locus — no need to merge disconnected HSPs heuristically.
    It is also substantially faster than tBLASTn on large genomes and directly
    outputs the spliced CDS sequence, so we skip the blastdbcmd extraction step.

QUERY PROTEIN:
    UniProt P69249 — rbcS from Nicotiana tabacum (180 aa, includes transit
    peptide). Downloaded at runtime from UniProt so it's always current.

WORKFLOW PER GENOME:
    1. Unzip the genome package downloaded by download_hifi_angiosperm_genomes.py.
    2. Locate the .fna genome FASTA inside the zip.
    3. Run miniprot with --gff to get a GFF3 + embedded spliced CDS output.
    4. Parse the GFF3: one locus per mRNA feature.
    5. Extract the spliced CDS sequence for each locus from the genome FASTA
       using the CDS coordinates (miniprot does not write a separate CDS fasta).
    6. Write per-genome hit FASTAs and a summary TSV across all genomes.

REQUIREMENTS:
    miniprot at /programs/miniprot-0.13/miniprot  (BioHPC path)
    Python standard library only (no third-party packages).

OUTPUT FILES (all in rbcS_results/):
    hits_summary.tsv          — one row per hit across all genomes
    <accession>_rbcS_hits.fna — spliced CDS + flanking context per genome
"""

import subprocess
import shutil
import zipfile
import os
import csv
import urllib.request
import glob

# ----------------------------- configuration ------------------------------- #

# Directory containing the genome zip files from the download script.
genome_zip_directory = "data/hifi_angiosperm_genomes"

# Where to write all results.
results_directory = "rbcS_results"

# UniProt accession for the rbcS query protein (Nicotiana tabacum rbcS).
# Ignored if reference_fasta_path is set below.
query_uniprot_accession = "P69249"

# Path to a local reference protein FASTA to use instead of downloading from
# UniProt. Set to None to download P69249 automatically.
# e.g. reference_fasta_path = "data/rca_reference.fa"
reference_fasta_path = "data/reference_protein_sequences/rca_reference.fa"

# Full path to miniprot binary (BioHPC).
miniprot_binary = "/programs/miniprot-0.13/miniprot"

# Number of CPU threads for miniprot.
miniprot_threads = 4

# Minimum fraction of the query protein that must be aligned to keep a hit.
# rbcS is 180 aa; 0.5 means at least 90 aa must be covered.
# Filters out tiny fragments while keeping partial copies.
minimum_query_coverage_fraction = 0.5

# Number of nucleotides of genomic sequence to extract on each side of the
# locus for flanking context — useful for inspecting upstream promoters etc.
flanking_nucleotides = 500

# --------------------------------------------------------------------------- #


def fetch_query_protein_fasta(uniprot_accession, output_fasta_path):
    """
    Download the query protein sequence from UniProt in FASTA format.
    Saves to output_fasta_path. Raises if the file is empty.
    """
    uniprot_fasta_url = (
        f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.fasta"
    )
    urllib.request.urlretrieve(uniprot_fasta_url, output_fasta_path)
    assert os.path.getsize(output_fasta_path) > 0, (
        f"Downloaded FASTA for {uniprot_accession} is empty."
    )
    print(f"  Query protein saved to: {output_fasta_path}")


def get_query_protein_length(query_fasta_path):
    """Return the number of amino acids in the first sequence of a FASTA file."""
    total_aa = 0
    with open(query_fasta_path) as fasta_file:
        for line in fasta_file:
            if not line.startswith(">"):
                total_aa += len(line.strip())
    assert total_aa > 0, f"Could not read sequence length from {query_fasta_path}"
    return total_aa


def find_genome_fasta_inside_zip(zip_path):
    """
    Return the zip-internal path of the largest *genomic.fna file inside an
    NCBI datasets genome zip. The largest file is used as a proxy for the
    primary assembly when multiple .fna files are present.
    """
    with zipfile.ZipFile(zip_path) as genome_zip:
        fna_entries = [
            entry for entry in genome_zip.namelist()
            if entry.endswith(".fna") and "genomic" in entry
        ]
    assert len(fna_entries) >= 1, (
        f"No *genomic.fna found in {zip_path}. "
        f"Contents: {zipfile.ZipFile(zip_path).namelist()[:10]}"
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


def run_miniprot(query_fasta_path, genome_fasta_path, gff_output_path, threads):
    """
    Run miniprot in GFF3 output mode.

    Key flags:
        --gff    : write GFF3 (includes ##PAF lines and CDS features)
        --outs   : minimum output score fraction; 0.9 keeps near-best hits
                   per locus (paralogs each get their own mRNA record)
        -t       : threads

    miniprot writes GFF3 to stdout; we capture it to gff_output_path.
    """
    with open(gff_output_path, "w") as gff_file:
        completed = subprocess.run(
            [
                miniprot_binary,
                "--gff",
                "--outs=0.9",   # report hits scoring ≥90% of the best hit
                "-t", str(threads),
                genome_fasta_path,
                query_fasta_path,
            ],
            stdout=gff_file,
            stderr=subprocess.PIPE,
            text=True,
        )
    assert completed.returncode == 0, (
        f"miniprot failed:\n{completed.stderr}"
    )


def parse_miniprot_gff3(gff_path, query_length_aa, assembly_accession):
    """
    Parse miniprot GFF3 output into a list of locus dicts.

    miniprot GFF3 structure:
        ##PAF  <tab-separated alignment summary line>
        <contig>  miniprot  mRNA  <start>  <end>  <score>  <strand>  .  ID=<id>;...
        <contig>  miniprot  CDS   <start>  <end>  .        <strand>  <phase>  Parent=<id>;...
        ...

    We collect one dict per mRNA with:
        locus_id, contig, genomic_start, genomic_end, strand,
        score, query_coverage_fraction,
        cds_intervals  (list of (start, end) in GFF3 1-based coords)

    Loci whose query coverage is below minimum_query_coverage_fraction are dropped.
    """
    loci = {}           # locus_id -> dict
    cds_by_locus = {}   # locus_id -> list of (start, end) intervals

    minimum_aligned_aa = query_length_aa * minimum_query_coverage_fraction

    with open(gff_path) as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue

            gff_columns = line.rstrip("\n").split("\t")
            if len(gff_columns) < 9:
                continue

            contig, source, feature_type, start_str, end_str, score_str, strand, phase, attributes = gff_columns
            start = int(start_str)
            end   = int(end_str)

            if feature_type == "mRNA":
                # Parse the attribute string into a dict.
                # GFF3 attributes look like: ID=MP000001;Identity=0.95;...
                attribute_dict = {}
                for attribute_pair in attributes.split(";"):
                    attribute_pair = attribute_pair.strip()
                    if "=" in attribute_pair:
                        key, value = attribute_pair.split("=", 1)
                        attribute_dict[key] = value

                locus_id = attribute_dict.get("ID", f"locus_{len(loci)}")

                # miniprot reports aligned length in the Identity/Rank fields;
                # the most direct coverage proxy is (end - start) vs query aa.
                # More accurately, we use the Target attribute if present:
                # Target=query <qstart> <qend>
                query_coverage_fraction = 0.0
                target_field = attribute_dict.get("Target", "")
                if target_field:
                    target_parts = target_field.split()
                    if len(target_parts) == 3:
                        aligned_aa = int(target_parts[2]) - int(target_parts[1])
                        query_coverage_fraction = aligned_aa / query_length_aa

                loci[locus_id] = {
                    "locus_id":                locus_id,
                    "contig":                  contig,
                    "genomic_start":           start,
                    "genomic_end":             end,
                    "strand":                  strand,
                    "score":                   score_str,
                    "query_coverage_fraction": round(query_coverage_fraction, 3),
                    "cds_intervals":           [],
                }

            elif feature_type == "CDS":
                # Find the parent mRNA ID.
                parent_id = None
                for attribute_pair in attributes.split(";"):
                    attribute_pair = attribute_pair.strip()
                    if attribute_pair.startswith("Parent="):
                        parent_id = attribute_pair.split("=", 1)[1]
                        break
                if parent_id and parent_id in loci:
                    loci[parent_id]["cds_intervals"].append((start, end))

    # Filter by query coverage and sort CDS intervals.
    passing_loci = []
    for locus in loci.values():
        if locus["query_coverage_fraction"] >= minimum_query_coverage_fraction:
            # Sort CDS intervals by position (important for correct CDS reconstruction).
            locus["cds_intervals"].sort(key=lambda interval: interval[0])
            passing_loci.append(locus)

    return passing_loci


def load_genome_sequences(genome_fasta_path):
    """
    Load a genome FASTA into a dict: {contig_id: sequence_string}.

    We load the whole genome into memory here. For chromosome-scale assemblies
    this can be a few GB. If memory is a concern, replace with an indexed
    approach (e.g. samtools faidx), but for the typical use case the simplicity
    is worth it.
    """
    sequences = {}
    current_id = None
    current_parts = []

    with open(genome_fasta_path) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_parts)
                # The contig ID is everything after ">" up to the first space.
                current_id = line[1:].split()[0]
                current_parts = []
            elif line:
                current_parts.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_parts)

    assert len(sequences) > 0, f"No sequences loaded from {genome_fasta_path}"
    return sequences


def extract_spliced_cds_with_flanking(locus, genome_sequences, assembly_accession,
                                       copy_number):
    """
    Reconstruct the spliced CDS nucleotide sequence by concatenating the CDS
    intervals from the genome, then wrap it with flanking genomic context.

    The flanking context is taken from the full genomic locus span (not per-exon),
    so the output sequence is:
        [flanking_upstream] + [full genomic locus incl. introns] + [flanking_downstream]

    This is what Tiberius expects — it needs to see the intron-containing locus
    in order to re-predict the gene model.

    Returns a tuple of (fasta_header, sequence_string) or None if the contig
    is not found in the genome.
    """
    contig_sequence = genome_sequences.get(locus["contig"])
    if contig_sequence is None:
        print(f"    WARNING: contig {locus['contig']} not found in genome FASTA.")
        return None

    contig_length = len(contig_sequence)

    # GFF3 is 1-based inclusive; Python slicing is 0-based exclusive.
    # Add flanking context, clamping to contig boundaries.
    extract_start_0based = max(0, locus["genomic_start"] - 1 - flanking_nucleotides)
    extract_end_0based   = min(contig_length, locus["genomic_end"] + flanking_nucleotides)

    extracted_sequence = contig_sequence[extract_start_0based:extract_end_0based]

    fasta_header = (
        f"{assembly_accession}_rbcS_copy{copy_number} "
        f"contig={locus['contig']} "
        f"locus={locus['genomic_start']}-{locus['genomic_end']} "
        f"strand={locus['strand']} "
        f"extracted={extract_start_0based + 1}-{extract_end_0based} "
        f"query_coverage={locus['query_coverage_fraction']:.2f} "
        f"n_exons={len(locus['cds_intervals'])}"
    )

    return fasta_header, extracted_sequence


def process_one_genome(zip_path, query_fasta_path, query_length_aa,
                       working_directory, results_directory):
    """
    Full miniprot pipeline for one genome zip. Returns summary rows (one per locus).
    """
    assembly_accession = os.path.basename(zip_path).replace(".zip", "")
    genome_work_dir = os.path.join(working_directory, assembly_accession)
    os.makedirs(genome_work_dir, exist_ok=True)

    print(f"\n  [{assembly_accession}] Extracting genome FASTA...")
    fna_entry = find_genome_fasta_inside_zip(zip_path)
    genome_fasta_path = extract_genome_fasta(zip_path, fna_entry, genome_work_dir)

    print(f"  [{assembly_accession}] Running miniprot...")
    gff_output_path = os.path.join(genome_work_dir, "miniprot_hits.gff")
    run_miniprot(query_fasta_path, genome_fasta_path, gff_output_path, miniprot_threads)

    loci = parse_miniprot_gff3(gff_output_path, query_length_aa, assembly_accession)
    print(f"  [{assembly_accession}] Loci passing coverage filter: {len(loci)}")

    if not loci:
        return []

    print(f"  [{assembly_accession}] Loading genome into memory for sequence extraction...")
    genome_sequences = load_genome_sequences(genome_fasta_path)

    hit_fasta_path = os.path.join(results_directory, f"{assembly_accession}_rbcS_hits.fna")
    summary_rows = []

    with open(hit_fasta_path, "w") as hit_fasta:
        for copy_number, locus in enumerate(loci, start=1):
            result = extract_spliced_cds_with_flanking(
                locus, genome_sequences, assembly_accession, copy_number
            )
            if result is None:
                continue
            fasta_header, sequence = result
            hit_fasta.write(f">{fasta_header}\n{sequence}\n")

            summary_rows.append({
                "assembly_accession":      assembly_accession,
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
    assert os.path.exists(miniprot_binary), (
        f"miniprot not found at {miniprot_binary}. "
        "Update miniprot_binary in the configuration section."
    )

    os.makedirs(results_directory, exist_ok=True)
    working_directory = os.path.join(results_directory, "working")
    os.makedirs(working_directory, exist_ok=True)

    # Step 1: get the query protein — either from a local file or UniProt.
    if reference_fasta_path is not None:
        assert os.path.exists(reference_fasta_path), (
            f"Reference FASTA not found: {reference_fasta_path}"
        )
        query_fasta_path = reference_fasta_path
        print(f"Using reference FASTA: {query_fasta_path}")
    else:
        query_fasta_path = os.path.join(results_directory, f"{query_uniprot_accession}.fasta")
        if not os.path.exists(query_fasta_path):
            print(f"Fetching query protein {query_uniprot_accession} from UniProt...")
            fetch_query_protein_fasta(query_uniprot_accession, query_fasta_path)
        else:
            print(f"Using cached query protein: {query_fasta_path}")

    query_length_aa = get_query_protein_length(query_fasta_path)
    print(f"Query protein length: {query_length_aa} aa")

    # Step 2: find all genome zips.
    genome_zip_paths = sorted(
        glob.glob(os.path.join(genome_zip_directory, "*.zip"))
    )
    assert len(genome_zip_paths) > 0, (
        f"No .zip files found in {genome_zip_directory}. "
        "Run download_hifi_angiosperm_genomes.py first."
    )
    print(f"Found {len(genome_zip_paths)} genome zip files to process.\n")

    # Step 3: process each genome and collect summary rows.
    all_summary_rows = []
    for zip_path in genome_zip_paths:
        print(f"Processing: {zip_path}")
        summary_rows = process_one_genome(
            zip_path, query_fasta_path, query_length_aa,
            working_directory, results_directory,
        )
        all_summary_rows.extend(summary_rows)

    # Step 4: write master summary TSV.
    summary_tsv_path = os.path.join(results_directory, "hits_summary.tsv")
    tsv_column_names = [
        "assembly_accession", "copy_number", "contig",
        "genomic_start", "genomic_end", "strand",
        "locus_length_nt", "n_exons",
        "query_coverage_fraction", "miniprot_score",
    ]
    with open(summary_tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_column_names, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_summary_rows)

    total_copies = len(all_summary_rows)
    genomes_with_hits = len({r["assembly_accession"] for r in all_summary_rows})
    print(f"\nDone.")
    print(f"  Genomes with at least one hit: {genomes_with_hits} / {len(genome_zip_paths)}")
    print(f"  Total rbcS loci found:         {total_copies}")
    print(f"  Summary TSV:  {summary_tsv_path}")
    print(f"  Per-genome FASTAs: {results_directory}/<accession>_rbcS_hits.fna")


if __name__ == "__main__":
    main()