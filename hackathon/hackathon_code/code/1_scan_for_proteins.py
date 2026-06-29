#!/usr/bin/env python3
"""
1_scan_for_hits.py

Step 1 of 3. Run on the genomics server (no GPU required).

Uses miniprot to scan angiosperm genome assemblies for copies of one or more
reference genes, extracts the genomic locus sequences and predicted peptides
for each hit, and bundles everything needed for folding into a single tar.gz
ready to copy to the GPU server.

OUTPUT (all under results_directory/):
    all_hits.fna              — genomic locus sequences (one per copy)
    all_hits_peps.faa         — predicted amino-acid sequences (one per copy)
    hits_summary.tsv          — one row per hit with coordinates and QC metadata
    for_folding.tar.gz        — all_hits_peps.faa + hits_summary.tsv, ready to
                                 scp to the GPU server for step 2

USAGE EXAMPLE:
    python 1_scan_for_hits.py
    scp results/gene_hits/for_folding.tar.gz gpu_server:~/rbcS_folding/
"""

import subprocess
import zipfile
import os
import csv
import glob
import tarfile


# Standard genetic code (codon → single-letter amino acid).
# '*' marks stop codons.
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(nucleotide_sequence):
    """Return the reverse complement of a DNA sequence."""
    return nucleotide_sequence.translate(COMPLEMENT)[::-1]


def translate_cds_intervals(cds_intervals, strand, contig_sequence):
    """
    Concatenate CDS exon intervals and translate to amino acids.

    GFF3 coordinates are 1-based inclusive. For a + strand gene, exons are
    concatenated in ascending genomic order. For a - strand gene, exons are
    concatenated in DESCENDING genomic order and then reverse-complemented as
    a unit before translation (the individual exons are on the minus strand).

    cds_intervals : list of (start, end) in 1-based inclusive GFF3 coords,
                    already sorted in ASCENDING genomic order by the parser.
    strand        : '+' or '-'
    contig_sequence : full contig nucleotide string (0-based Python indexing)

    Returns the amino acid string, stopping at (and excluding) the first stop
    codon. Returns an empty string if translation fails for any reason.
    """
    # Extract each exon as a nucleotide string (convert to 0-based slice).
    exon_sequences = []
    for (gff_start, gff_end) in cds_intervals:
        # GFF3: start and end are 1-based inclusive.
        # Python slice: subtract 1 from start; end is already exclusive.
        exon_sequences.append(contig_sequence[gff_start - 1 : gff_end])

    if strand == '-':
        # Reverse the exon order and reverse-complement the concatenated sequence.
        # (Each exon is on the minus strand; concatenating in reverse genomic order
        # and taking the RC gives the correct coding sequence 5'→3'.)
        concatenated_cds = reverse_complement("".join(reversed(exon_sequences)))
    else:
        concatenated_cds = "".join(exon_sequences)

    concatenated_cds = concatenated_cds.upper()

    # Translate codon by codon; stop at the first stop codon.
    amino_acids = []
    for codon_start in range(0, len(concatenated_cds) - 2, 3):
        codon = concatenated_cds[codon_start : codon_start + 3]
        amino_acid = CODON_TABLE.get(codon, 'X')   # 'X' for any ambiguous codon
        if amino_acid == '*':
            break
        amino_acids.append(amino_acid)

    return "".join(amino_acids)


# ----------------------------- configuration ------------------------------- #

# Genome zip files downloaded from NCBI datasets.
genome_zip_directory = "data/hifi_angiosperm_genomes"

# One or more reference protein FASTAs (.fa / .faa / .fasta).
# All files here are concatenated into a single miniprot query.
reference_proteins_directory = "data/reference_protein_sequences"

# Where to write results.
results_directory = "results/gene_hits_scan"

# Full path to miniprot binary.
miniprot_binary = "/programs/miniprot-0.13/miniprot"

# Number of CPU threads for miniprot.
miniprot_threads = 10

# Minimum fraction of the reference protein that must be aligned to keep a hit.
# 0.5 = at least half the reference must be covered.
minimum_query_coverage_fraction = 0.5

# Nucleotides of genomic context on each side of each hit locus (for Tiberius).
flanking_nucleotides = 1500

# --------------------------------------------------------------------------- #


def collect_and_concatenate_reference_proteins(reference_proteins_directory,
                                               output_fasta_path):
    """
    Concatenate all reference FASTAs into one combined query file.
    Returns {protein_id: length_in_aa} used later to compute query coverage.
    """
    fasta_extensions = ("*.fa", "*.faa", "*.fasta")
    reference_fasta_paths = []
    for extension in fasta_extensions:
        reference_fasta_paths.extend(
            glob.glob(os.path.join(reference_proteins_directory, extension))
        )
    reference_fasta_paths = sorted(reference_fasta_paths)

    assert len(reference_fasta_paths) > 0, (
        f"No FASTA files found in {reference_proteins_directory}."
    )
    print(f"Found {len(reference_fasta_paths)} reference protein file(s):")
    for path in reference_fasta_paths:
        print(f"  {path}")

    protein_lengths_by_id = {}
    with open(output_fasta_path, "w") as combined_fasta:
        for fasta_path in reference_fasta_paths:
            current_protein_id = None
            current_length = 0
            with open(fasta_path) as individual_fasta:
                for line in individual_fasta:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if current_protein_id is not None:
                            protein_lengths_by_id[current_protein_id] = current_length
                        current_protein_id = line[1:].split()[0]
                        current_length = 0
                    else:
                        current_length += len(line)
                    combined_fasta.write(line + "\n")
                if current_protein_id is not None:
                    protein_lengths_by_id[current_protein_id] = current_length

    assert len(protein_lengths_by_id) > 0, "No protein sequences were read."
    print(f"Total reference proteins: {len(protein_lengths_by_id)}")
    return protein_lengths_by_id


def find_genome_fasta_inside_zip(zip_path):
    """Return the zip-internal path of the largest *genomic.fna inside an NCBI datasets zip."""
    with zipfile.ZipFile(zip_path) as genome_zip:
        fna_entries = [
            entry for entry in genome_zip.namelist()
            if entry.endswith(".fna") and "genomic" in entry
        ]
    assert len(fna_entries) >= 1, f"No *genomic.fna found in {zip_path}."
    if len(fna_entries) > 1:
        with zipfile.ZipFile(zip_path) as genome_zip:
            fna_entries = sorted(
                fna_entries,
                key=lambda entry: genome_zip.getinfo(entry).file_size,
                reverse=True,
            )
    return fna_entries[0]


def extract_genome_fasta(zip_path, fna_entry_name, extract_to_directory):
    """Extract one .fna entry from the zip and return the path of the extracted file."""
    with zipfile.ZipFile(zip_path) as genome_zip:
        genome_zip.extract(fna_entry_name, extract_to_directory)
    extracted_path = os.path.join(extract_to_directory, fna_entry_name)
    assert os.path.exists(extracted_path), f"Extraction failed: {extracted_path}"
    return extracted_path


def run_miniprot(combined_query_fasta_path, genome_fasta_path,
                 gff_output_path, threads):
    """
    Run miniprot with all reference proteins against one genome.

    --gff   : GFF3 output (one mRNA + CDS features per hit locus)
    --outs  : retain hits scoring >= 90% of the top hit per locus so paralogs
              are not suppressed by a best-hit filter
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
    Parse miniprot GFF3 into a list of hit-locus dicts.

    The Target= attribute on each mRNA line records which reference protein
    matched and the aligned amino-acid range:
        Target=<protein_id> <aa_start> <aa_end>

    query_coverage_fraction = aligned_aa / total_reference_protein_length
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

            attribute_dict = {}
            for pair in attributes.split(";"):
                pair = pair.strip()
                if "=" in pair:
                    key, value = pair.split("=", 1)
                    attribute_dict[key] = value

            if feature_type == "mRNA":
                locus_id = attribute_dict.get("ID", f"locus_{len(loci_by_id)}")
                gene_name = "unknown"
                query_coverage_fraction = 0.0
                target_field = attribute_dict.get("Target", "")
                if target_field:
                    target_parts = target_field.split()
                    if len(target_parts) == 3:
                        gene_name = target_parts[0]
                        aligned_aa = int(target_parts[2]) - int(target_parts[1])
                        protein_length = protein_lengths_by_id.get(gene_name, 1)
                        query_coverage_fraction = aligned_aa / protein_length

                loci_by_id[locus_id] = {
                    "locus_id": locus_id,
                    "gene_name": gene_name,
                    "contig": contig,
                    "genomic_start": int(start_str),
                    "genomic_end": int(end_str),
                    "strand": strand,
                    "score": score_str,
                    "query_coverage_fraction": round(query_coverage_fraction, 3),
                    "cds_intervals": [],
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
    """Read a genome FASTA into {contig_id: nucleotide_sequence}."""
    sequences = {}
    current_id = None
    current_parts = []
    with open(genome_fasta_path) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_parts)
                current_id = line[1:].split()[0]
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
    Extract the full genomic span of one hit locus plus flanking_nucleotides
    of context on each side.

    GFF3 coordinates are 1-based inclusive; Python slicing is 0-based exclusive,
    so we subtract 1 from the start and clamp both ends to the contig length.

    Header format: {accession}_{gene_name}_copy{N}
    Returns (fasta_header, sequence) or None if the contig is missing.
    """
    contig_sequence = genome_sequences.get(locus["contig"])
    if contig_sequence is None:
        print(f"    WARNING: contig {locus['contig']} not found — skipping.")
        return None

    contig_length = len(contig_sequence)
    extract_start = max(0, locus["genomic_start"] - 1 - flanking_nucleotides)
    extract_end = min(contig_length, locus["genomic_end"] + flanking_nucleotides)
    extracted_sequence = contig_sequence[extract_start:extract_end]

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
                       working_directory, all_hits_fna_handle, all_hits_peps_handle):
    """
    Full pipeline for one genome zip:
      extract → miniprot → parse GFF → extract sequences → write to shared FASTAs.

    Sequences are appended to open file handles so only one genome is in memory
    at a time. Returns a list of summary row dicts (one per passing hit).
    """
    assembly_accession = os.path.basename(zip_path).replace(".zip", "")
    genome_work_dir = os.path.join(working_directory, assembly_accession)
    os.makedirs(genome_work_dir, exist_ok=True)

    print(f"\n  [{assembly_accession}] Extracting genome FASTA...")
    fna_entry = find_genome_fasta_inside_zip(zip_path)
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

    copy_counter_by_gene = {}
    summary_rows = []

    for locus in passing_loci:
        gene_name = locus["gene_name"]
        copy_number = copy_counter_by_gene.get(gene_name, 0) + 1
        copy_counter_by_gene[gene_name] = copy_number

        # Unique copy ID shared between the nucleotide and peptide FASTAs.
        safe_gene_name = gene_name.replace("|", "_").replace(" ", "_")
        copy_id = f"{assembly_accession}_{safe_gene_name}_copy{copy_number}"

        # --- nucleotide locus ---
        result = extract_genomic_locus_with_flanking(
            locus, genome_sequences, assembly_accession, copy_number
        )
        if result is None:
            continue
        fasta_header, locus_sequence = result
        all_hits_fna_handle.write(f">{fasta_header}\n{locus_sequence}\n")

        # --- predicted peptide (translated from CDS intervals) ---
        # Translate the CDS exons directly from the genome sequence.
        # This is necessary because miniprot does not emit an MP= tag in GFF3
        # output unless --aln is passed; the CDS coordinates are always present.
        peptide_sequence = translate_cds_intervals(
            locus["cds_intervals"], locus["strand"],
            genome_sequences[locus["contig"]]
        )
        peptide_header = (
            f"{copy_id} "
            f"gene={gene_name} "
            f"contig={locus['contig']} "
            f"locus={locus['genomic_start']}-{locus['genomic_end']} "
            f"strand={locus['strand']} "
            f"query_coverage={locus['query_coverage_fraction']:.2f}"
        )
        if peptide_sequence:
            all_hits_peps_handle.write(f">{peptide_header}\n{peptide_sequence}\n")
        else:
            print(f"    WARNING: translation produced empty peptide for {copy_id}.")

        summary_rows.append({
            "copy_id": copy_id,
            "assembly_accession": assembly_accession,
            "gene_name": gene_name,
            "copy_number": copy_number,
            "contig": locus["contig"],
            "genomic_start": locus["genomic_start"],
            "genomic_end": locus["genomic_end"],
            "strand": locus["strand"],
            "locus_length_nt": locus["genomic_end"] - locus["genomic_start"] + 1,
            "n_exons": len(locus["cds_intervals"]),
            "query_coverage_fraction": locus["query_coverage_fraction"],
            "miniprot_score": locus["score"],
            "has_predicted_peptide": bool(peptide_sequence),
        })

    return summary_rows


def bundle_for_transit(results_directory):
    """
    Create for_folding.tar.gz containing just the files needed on the GPU server:
        all_hits_peps.faa  — peptides to fold
        hits_summary.tsv   — metadata (copy_id links results back to this table)
    """
    bundle_path = os.path.join(results_directory, "for_folding.tar.gz")
    files_to_bundle = [
        os.path.join(results_directory, "all_hits_peps.faa"),
        os.path.join(results_directory, "hits_summary.tsv"),
    ]
    for file_path in files_to_bundle:
        assert os.path.exists(file_path), f"Expected output file missing: {file_path}"

    with tarfile.open(bundle_path, "w:gz") as tar:
        for file_path in files_to_bundle:
            # Store with basename only so extraction is flat.
            tar.add(file_path, arcname=os.path.basename(file_path))

    bundle_size_mb = os.path.getsize(bundle_path) / 1e6
    print(f"\n[bundle] {bundle_path}  ({bundle_size_mb:.1f} MB)")
    print("[bundle] copy to GPU server with:")
    print(f"    scp {bundle_path} gpu_server:~/rbcS_folding/")
    return bundle_path


def main():
    assert os.path.exists(miniprot_binary), \
        f"miniprot not found at {miniprot_binary}."
    assert os.path.isdir(reference_proteins_directory), \
        f"Reference proteins directory not found: {reference_proteins_directory}"

    os.makedirs(results_directory, exist_ok=True)
    working_directory = os.path.join(results_directory, "working")
    os.makedirs(working_directory, exist_ok=True)

    combined_query_fasta_path = os.path.join(working_directory, "all_reference_proteins.faa")
    protein_lengths_by_id = collect_and_concatenate_reference_proteins(
        reference_proteins_directory, combined_query_fasta_path
    )

    genome_zip_paths = sorted(glob.glob(os.path.join(genome_zip_directory, "*.zip")))
    assert len(genome_zip_paths) > 0, \
        f"No .zip files found in {genome_zip_directory}."
    print(f"\nFound {len(genome_zip_paths)} genome zip files to process.")

    all_hits_fna_path = os.path.join(results_directory, "all_hits.fna")
    all_hits_peps_path = os.path.join(results_directory, "all_hits_peps.faa")
    summary_tsv_path = os.path.join(results_directory, "hits_summary.tsv")

    all_summary_rows = []

    with open(all_hits_fna_path, "w") as all_hits_fna_handle, \
         open(all_hits_peps_path, "w") as all_hits_peps_handle:
        for zip_path in genome_zip_paths:
            print(f"\nProcessing: {zip_path}")
            summary_rows = process_one_genome(
                zip_path, combined_query_fasta_path, protein_lengths_by_id,
                working_directory, all_hits_fna_handle, all_hits_peps_handle,
            )
            all_summary_rows.extend(summary_rows)

    tsv_column_names = [
        "copy_id", "assembly_accession", "gene_name", "copy_number", "contig",
        "genomic_start", "genomic_end", "strand",
        "locus_length_nt", "n_exons",
        "query_coverage_fraction", "miniprot_score", "has_predicted_peptide",
    ]
    with open(summary_tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_column_names, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_summary_rows)

    total_hits = len(all_summary_rows)
    genomes_with_hits = len({r["assembly_accession"] for r in all_summary_rows})
    genes_found = sorted({r["gene_name"] for r in all_summary_rows})
    copies_with_peptides = sum(r["has_predicted_peptide"] for r in all_summary_rows)

    print(f"\nDone.")
    print(f"  Genomes with at least one hit : {genomes_with_hits} / {len(genome_zip_paths)}")
    print(f"  Total hit loci                : {total_hits}")
    print(f"  Copies with predicted peptide : {copies_with_peptides}")
    print(f"  Genes recovered               : {genes_found}")

    bundle_for_transit(results_directory)


if __name__ == "__main__":
    main()