#!/usr/bin/env python3
"""
3_build_qc_table.py

Step 3 of 3. Run on the genomics server after receiving folded_structures.tar.gz.

Unpacks the folded structures bundle, aligns every candidate CIF against a
single reference structure with US-align (reference-based, not all-vs-all),
aligns every candidate peptide against a single reference protein (Needleman-
Wunsch / BLOSUM62), and assembles one row per copy into a single QC table.

Nothing is dropped — every copy gets pass/fail flags so you can embed the
full set (ESM, PlantCAD2) and watch the QC filters stratify in embedding space.

USAGE:
    # Receive bundle from GPU server:
    scp gpu_server:~/rbcS_folding/folded_structures.tar.gz results/gene_hits/

    # Then run:
    python 3_build_qc_table.py

OUTPUT:
    results/rbcs_qc/rbcs_qc_table.tsv
"""

import csv
import math
import os
import re
import subprocess
import tarfile
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


# ----------------------------- configuration ------------------------------- #

# Results from step 1 (on this server).
scan_results_directory = "results/gene_hits_rbcs"
candidate_peptides_fasta_path = "results/gene_hits_rbcs/all_hits_peps.faa"
candidate_nucleotides_fasta_path = "results/gene_hits_rbcs/all_hits.fna"
hits_summary_tsv_path = "results/gene_hits_rbcs/hits_summary.tsv"

# Bundle received back from the GPU server.
folded_structures_bundle_path = "results/gene_hits_rbcs/folded_structures.tar.gz"

# Where the bundle will be unpacked.
folded_structures_directory = "results/gene_hits_rbcs/folded_structures"

# Reference protein (single sequence FASTA) for sequence-level alignment QC.
reference_peptide_fasta_path = "data/temp_ref_dir/1RCX_reference_proteins.faa"

# Reference structure for structural QC (CIF or PDB).
# If left as None, the script will look for REFERENCE.cif in the unpacked bundle
# (i.e. if you folded the reference on the GPU server alongside the candidates).
reference_structure_path = "data/structures/1RCX.cif"  # e.g. "data/reference_structures/rbcS_reference.cif"

# Where to write the final QC table.
output_directory = "results/rbcs_qc"

# Full path to US-align binary.
usalign_binary = "/home/hdd29/USalign/USalign"

# Chain to extract from the reference structure for structural QC.
# The rbcS crystal structure 1RCX is a hexadecamer (chains A-P); aligning a
# single-chain ESMFold2 model against all 16 chains inflates the normalisation
# length ~16-fold, producing near-zero TM-scores. Set this to any single rbcS
# subunit chain (e.g. "A") to get a meaningful per-subunit comparison.
# Set to None to use the reference structure as-is (correct for monomeric refs).
reference_chain_id = "P"

# QC thresholds.
# pid_threshold     : minimum percent identity to reference protein (0–1).
#                     0.30 is permissive for angiosperm rbcS; raise if needed.
# coverage_threshold: minimum fraction of BOTH query and reference that must be
#                     covered by the alignment (0–1). 0.5 = at least half covered.
# tm_threshold      : TM-score threshold for "same fold" (0–1). 0.5 is the
#                     conventional cutoff; normalised by reference length so
#                     truncated models score low even if their fragment is correct.
pid_threshold = 0.30
coverage_threshold = 0.50
tm_threshold = 0.50

# --------------------------------------------------------------------------- #

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def remove_nonstandard_residues(raw_sequence):
    """Uppercase, strip stop characters (*), keep only the 20 standard amino acids."""
    raw_sequence = str(raw_sequence).upper().replace("*", "")
    return "".join(residue for residue in raw_sequence if residue in STANDARD_AMINO_ACIDS)


def make_filesystem_safe_id(sequence_id):
    """Replace characters unsafe for filenames with underscores."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", sequence_id)


def read_fasta_as_dict(fasta_path):
    """Return {first_token_of_header: raw_sequence_string} for a FASTA file."""
    sequences_by_id = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences_by_id[record.id] = str(record.seq)
    assert len(sequences_by_id) > 0, f"No sequences read from {fasta_path}"
    return sequences_by_id


def unpack_folded_structures_bundle(bundle_path, extract_to_directory):
    """
    Unpack folded_structures.tar.gz into extract_to_directory.
    Expected contents: fold_metrics.tsv and cif/*.cif
    Returns the path to fold_metrics.tsv and the cif subdirectory.
    """
    assert os.path.exists(bundle_path), (
        f"Folded structures bundle not found: {bundle_path}\n"
        "Did you scp folded_structures.tar.gz back from the GPU server?"
    )
    os.makedirs(extract_to_directory, exist_ok=True)
    with tarfile.open(bundle_path, "r:gz") as tar:
        tar.extractall(extract_to_directory)

    metrics_tsv_path = Path(extract_to_directory) / "fold_metrics.tsv"
    cif_directory = Path(extract_to_directory) / "cif"
    assert metrics_tsv_path.exists(), f"fold_metrics.tsv missing from bundle."
    assert cif_directory.is_dir(), f"cif/ directory missing from bundle."
    return metrics_tsv_path, cif_directory


def read_fold_metrics(metrics_tsv_path):
    """
    Read the fold_metrics.tsv from step 2 into a dict keyed by copy_id.
    Casts numeric fields to float (NaN for missing values).
    """
    fold_metrics_by_copy_id = {}
    with open(metrics_tsv_path) as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            for numeric_column in ("plddt_mean", "ptm", "iptm"):
                try:
                    row[numeric_column] = float(row[numeric_column])
                except (ValueError, KeyError):
                    row[numeric_column] = math.nan
            row["fold_ok"] = row.get("fold_ok", "False") == "True"
            fold_metrics_by_copy_id[row["copy_id"]] = row
    return fold_metrics_by_copy_id


# ---------------------------------------------------------------------------
# Sequence-level QC: BLOSUM62 global alignment vs one reference protein
# ---------------------------------------------------------------------------

def build_protein_aligner():
    """
    Global aligner (Needleman-Wunsch) with BLOSUM62 and free end gaps.

    Free end gaps mean that terminal truncation of a candidate is not penalised
    in the score — important because many candidates are fragments. Coverage
    metrics still capture truncation explicitly.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11.0
    aligner.extend_gap_score = -1.0
    # Free end gaps: terminal gaps on either side do not contribute to the score.
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    return aligner


def compute_identity_and_coverage(reference_sequence, candidate_sequence, aligner):
    """
    Align candidate vs reference and return (percent_identity, coverage_self, coverage_ref).

    Definitions:
        n_aligned_pairs  = alignment columns where BOTH sequences have a residue
        n_identical      = aligned columns where both residues are the same

        percent_identity = n_identical    / n_aligned_pairs
        coverage_self    = n_aligned_pairs / len(candidate)   [fraction of candidate covered]
        coverage_ref     = n_aligned_pairs / len(reference)   [fraction of reference covered]

    All three are fractions in [0, 1]. Returns NaN for empty sequences.
    """
    if len(candidate_sequence) == 0:
        return math.nan, math.nan, math.nan

    best_alignment = aligner.align(reference_sequence, candidate_sequence)[0]

    # .aligned gives (reference_blocks, candidate_blocks) as lists of [start, end) pairs.
    # Each reference_blocks[i] aligns to candidate_blocks[i] with the same block length.
    reference_blocks, candidate_blocks = best_alignment.aligned

    number_of_aligned_pairs = 0
    number_of_identical_residues = 0
    for (ref_start, ref_end), (cand_start, cand_end) in zip(reference_blocks, candidate_blocks):
        block_length = ref_end - ref_start
        number_of_aligned_pairs += block_length
        for offset in range(block_length):
            if reference_sequence[ref_start + offset] == candidate_sequence[cand_start + offset]:
                number_of_identical_residues += 1

    if number_of_aligned_pairs == 0:
        return 0.0, 0.0, 0.0

    percent_identity = number_of_identical_residues / number_of_aligned_pairs
    coverage_self = number_of_aligned_pairs / len(candidate_sequence)
    coverage_ref = number_of_aligned_pairs / len(reference_sequence)
    return percent_identity, coverage_self, coverage_ref


# ---------------------------------------------------------------------------
# Structure-level QC: US-align vs one reference structure
# ---------------------------------------------------------------------------

def extract_single_chain_from_structure(input_structure_path, chain_id, output_path):
    """
    Write a new structure file containing only the specified chain.

    This handles the common case where the reference is a multi-chain crystal
    structure (e.g. the rbcS hexadecamer 1RCX has chains A-P). US-align aligns
    against ALL chains by default, which inflates the normalisation length and
    produces near-zero TM-scores for a single-chain candidate model.

    Works for both PDB (.pdb) and mmCIF (.cif) format by line-level filtering:
      - PDB : keep ATOM/HETATM lines where column 22 (0-based) == chain_id,
              and keep all REMARK/HEADER lines.
      - CIF : keep _atom_site lines where the auth_asym_id field == chain_id.
              For CIF this uses a simple heuristic: if the line has enough
              whitespace-delimited fields and the chain field matches, keep it.
              Also keeps all non-atom-site lines (header, cell parameters, etc.)
              so the file remains valid.

    If extraction yields zero ATOM records the function raises AssertionError.
    """
    input_path = Path(input_structure_path)
    output_path = Path(output_path)
    suffix = input_path.suffix.lower()

    atom_lines_written = 0

    if suffix == ".pdb":
        with open(input_path) as infile, open(output_path, "w") as outfile:
            for line in infile:
                record_type = line[:6].strip()
                if record_type in ("ATOM", "HETATM"):
                    # PDB format: chain ID is character at column index 21 (1-based col 22).
                    if len(line) > 21 and line[21] == chain_id:
                        outfile.write(line)
                        atom_lines_written += 1
                else:
                    # Keep all non-coordinate lines (HEADER, REMARK, END, etc.)
                    outfile.write(line)

    elif suffix == ".cif":
        # mmCIF: we need to identify which column index holds auth_asym_id.
        # Parse the loop_ block column headers to find it dynamically.
        with open(input_path) as infile:
            lines = infile.readlines()

        in_atom_site_loop = False
        column_headers = []
        auth_asym_id_column_index = None
        output_lines = []

        for line in lines:
            stripped = line.strip()

            if stripped == "loop_":
                # Reset loop state; we may be entering a new loop block.
                in_atom_site_loop = False
                column_headers = []
                auth_asym_id_column_index = None
                output_lines.append(line)
                continue

            if stripped.startswith("_atom_site."):
                in_atom_site_loop = True
                column_headers.append(stripped)
                if stripped == "_atom_site.auth_asym_id":
                    auth_asym_id_column_index = len(column_headers) - 1
                output_lines.append(line)
                continue

            if in_atom_site_loop and auth_asym_id_column_index is not None:
                # This is a data row in the _atom_site loop.
                if stripped.startswith("_") or stripped in ("#", "loop_", ""):
                    # Leaving the atom site block.
                    in_atom_site_loop = False
                    output_lines.append(line)
                    continue
                fields = stripped.split()
                if len(fields) > auth_asym_id_column_index:
                    row_chain = fields[auth_asym_id_column_index]
                    if row_chain == chain_id:
                        output_lines.append(line)
                        atom_lines_written += 1
                    # Silently skip rows for other chains.
                else:
                    output_lines.append(line)
                continue

            output_lines.append(line)

        with open(output_path, "w") as outfile:
            outfile.writelines(output_lines)

    else:
        raise ValueError(f"Unsupported structure format: {suffix} (expected .pdb or .cif)")

    assert atom_lines_written > 0, (
        f"No ATOM records for chain {chain_id} found in {input_structure_path}. "
        "Check that reference_chain_id matches a real chain in the structure."
    )
    return output_path


def parse_usalign_output(usalign_text):
    """
    Extract both TM-scores from US-align stdout.

    US-align prints TM-scores in this order:
        TM-score=X  (normalised by length of structure 1, i.e. the candidate)
        TM-score=Y  (normalised by length of structure 2, i.e. the reference)

    We call USalign <candidate> <reference>, so the second score asks:
    "What fraction of the reference fold does this candidate reproduce?"
    """
    tm_scores = [float(v) for v in re.findall(r"TM-score=\s*([0-9.]+)", usalign_text)]
    tm_normalised_by_candidate = tm_scores[0] if len(tm_scores) > 0 else math.nan
    tm_normalised_by_reference = tm_scores[1] if len(tm_scores) > 1 else math.nan
    return tm_normalised_by_candidate, tm_normalised_by_reference


def compute_tm_score_vs_reference(candidate_cif_path, reference_single_chain_path):
    """
    Run US-align (candidate first, reference second) and return the TM-score
    normalised by reference length. NaN on failure.

    reference_single_chain_path should already be a single-chain structure
    (produced by extract_single_chain_from_structure) so the normalisation
    length matches one rbcS subunit, not the whole hexadecamer.
    """
    superposed_output_path = candidate_cif_path.parent / f"{candidate_cif_path.stem}_superposed"

    completed = subprocess.run(
        [usalign_binary, str(candidate_cif_path), str(reference_single_chain_path),
        "-o", str(superposed_output_path)],
        capture_output=True, text=True,
    )
    if completed.returncode != 0:
        print(f"    WARNING: US-align failed for {candidate_cif_path.name}: {completed.stderr[:200]}")
        return math.nan
    _, tm_normalised_by_reference = parse_usalign_output(
        completed.stdout + "\n" + completed.stderr
    )
    return tm_normalised_by_reference


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # --- 0. validate inputs -------------------------------------------------
    for required_path in [
        candidate_peptides_fasta_path,
        candidate_nucleotides_fasta_path,
        reference_peptide_fasta_path,
        folded_structures_bundle_path,
    ]:
        assert os.path.exists(required_path), f"Required input not found: {required_path}"

    os.makedirs(output_directory, exist_ok=True)

    # --- 1. unpack folded structures ----------------------------------------
    print("[1] unpacking folded structures bundle...")
    fold_metrics_tsv_path, cif_directory = unpack_folded_structures_bundle(
        folded_structures_bundle_path, folded_structures_directory
    )
    fold_metrics_by_copy_id = read_fold_metrics(fold_metrics_tsv_path)
    print(f"    {len(fold_metrics_by_copy_id)} fold metric rows loaded")

    # --- 2. locate reference structure and extract single chain if needed ----
    resolved_reference_structure_path = reference_structure_path
    if resolved_reference_structure_path is None:
        # Fall back to a REFERENCE.cif that may have been folded on the GPU server.
        candidate_reference_cif = cif_directory / "REFERENCE.cif"
        assert candidate_reference_cif.exists(), (
            "No reference_structure_path set and REFERENCE.cif not found in bundle.\n"
            "Either set reference_structure_path at the top of this script, or\n"
            "re-run step 2 with the reference peptide included in all_hits_peps.faa."
        )
        resolved_reference_structure_path = str(candidate_reference_cif)

    # If the reference is a multi-chain assembly (e.g. a hexadecamer crystal
    # structure), extract just one chain so US-align normalises by a single
    # subunit length rather than the whole complex.
    if reference_chain_id is not None:
        reference_single_chain_path = (
            Path(output_directory) /
            f"reference_chain_{reference_chain_id}{Path(resolved_reference_structure_path).suffix}"
        )
        if not reference_single_chain_path.exists():
            print(f"[2] extracting chain {reference_chain_id} from reference structure...")
            extract_single_chain_from_structure(
                resolved_reference_structure_path,
                reference_chain_id,
                reference_single_chain_path,
            )
        else:
            print(f"[2] reusing existing chain {reference_chain_id} extract: {reference_single_chain_path}")
    else:
        reference_single_chain_path = Path(resolved_reference_structure_path)

    print(f"[2] reference for US-align: {reference_single_chain_path}")

    # --- 3. load sequences --------------------------------------------------
    print("[3] loading candidate sequences...")
    candidate_peptides_by_id = read_fasta_as_dict(candidate_peptides_fasta_path)
    candidate_nucleotides_by_id = read_fasta_as_dict(candidate_nucleotides_fasta_path)
    print(f"    {len(candidate_peptides_by_id)} candidate peptides")

    reference_records = list(SeqIO.parse(reference_peptide_fasta_path, "fasta"))
    assert len(reference_records) >= 1, "Reference peptide FASTA is empty."
    reference_amino_acid_sequence = remove_nonstandard_residues(reference_records[0].seq)
    assert len(reference_amino_acid_sequence) > 0, "Reference peptide has no standard residues."
    print(f"    reference length: {len(reference_amino_acid_sequence)} aa")

    # --- 4. build the QC table ---------------------------------------------
    print("[4] computing per-copy QC metrics...")
    protein_aligner = build_protein_aligner()
    table_rows = []

    for copy_id in sorted(candidate_peptides_by_id):
        safe_id = make_filesystem_safe_id(copy_id)

        # Amino acid sequence.
        candidate_amino_acid_sequence = remove_nonstandard_residues(
            candidate_peptides_by_id[copy_id]
        )

        # Nucleotide sequence (matched by the same first-token ID).
        candidate_nucleotide_sequence = candidate_nucleotides_by_id.get(copy_id, "")
        if candidate_nucleotide_sequence == "":
            print(f"    WARNING: no nucleotide sequence found for {copy_id}")

        candidate_length = len(candidate_amino_acid_sequence)

        # --- sequence-level QC ---
        # Align candidate vs reference protein; compute identity and coverage.
        percent_identity, coverage_self, coverage_ref = compute_identity_and_coverage(
            reference_amino_acid_sequence, candidate_amino_acid_sequence, protein_aligner
        )
        # Pass if identity is sufficient AND the alignment covers at least
        # coverage_threshold of BOTH sequences (short fragments fail even if
        # the matched portion is identical).
        sequence_filter_passed = (
            candidate_length > 0
            and percent_identity >= pid_threshold
            and coverage_self >= coverage_threshold
            and coverage_ref >= coverage_threshold
        )

        # --- structure-level QC ---
        # Retrieve fold metrics from step 2 and run US-align vs reference.
        fold_row = fold_metrics_by_copy_id.get(copy_id, {})
        fold_ok = fold_row.get("fold_ok", False)
        plddt_mean = fold_row.get("plddt_mean", math.nan)
        ptm = fold_row.get("ptm", math.nan)
        iptm = fold_row.get("iptm", math.nan)

        tm_score = math.nan
        if fold_ok:
            candidate_cif_path = cif_directory / f"{safe_id}.cif"
            if candidate_cif_path.exists():
                tm_score = compute_tm_score_vs_reference(
                    candidate_cif_path,
                    reference_single_chain_path,
                )
            else:
                print(f"    WARNING: fold_ok=True but CIF not found for {copy_id}")

        # Pass the structure filter if TM-score meets the "same fold" threshold.
        structure_filter_passed = (not math.isnan(tm_score)) and (tm_score >= tm_threshold)

        table_rows.append({
            "species": copy_id.rsplit("_copy", 1)[0],   # accession portion of the ID
            "gene_id": copy_id,
            "nucleotide_seq": candidate_nucleotide_sequence,
            "amino_acid_seq": candidate_amino_acid_sequence,
            "pid_vs_ref": round(percent_identity, 4) if not math.isnan(percent_identity) else math.nan,
            "coverage_self": round(coverage_self, 4) if not math.isnan(coverage_self) else math.nan,
            "coverage_ref": round(coverage_ref, 4) if not math.isnan(coverage_ref) else math.nan,
            "plddt_mean": round(plddt_mean, 3) if not math.isnan(plddt_mean) else math.nan,
            "ptm": round(ptm, 4) if not math.isnan(ptm) else math.nan,
            "iptm": round(iptm, 4) if not math.isnan(iptm) else math.nan,
            "tm_score": round(tm_score, 4) if not math.isnan(tm_score) else math.nan,
            "length": candidate_length,
            "aln_filter_passed": sequence_filter_passed,
            "fold_filter_passed": structure_filter_passed,
        })

    # Sanity check: one row per input sequence (nothing dropped).
    assert len(table_rows) == len(candidate_peptides_by_id), \
        "Row count mismatch — a copy was lost during table construction."

    # --- 5. write table ----------------------------------------------------
    output_table_path = Path(output_directory) / "rbcs_qc_table.tsv"
    column_order = [
        "species", "gene_id",
        "length",
        "pid_vs_ref", "coverage_self", "coverage_ref",
        "plddt_mean", "ptm", "iptm",
        "tm_score",
        "aln_filter_passed", "fold_filter_passed",
        "amino_acid_seq", "nucleotide_seq",
    ]
    with open(output_table_path, "w") as table_file:
        table_file.write("\t".join(column_order) + "\n")
        for row in table_rows:
            table_file.write("\t".join(str(row[col]) for col in column_order) + "\n")

    number_passing_sequence = sum(r["aln_filter_passed"] for r in table_rows)
    number_passing_structure = sum(r["fold_filter_passed"] for r in table_rows)
    number_passing_both = sum(
        r["aln_filter_passed"] and r["fold_filter_passed"] for r in table_rows
    )
    print(f"\n[done] total copies             : {len(table_rows)}")
    print(f"[done] pass sequence filter     : {number_passing_sequence}")
    print(f"[done] pass structure filter    : {number_passing_structure}")
    print(f"[done] pass both filters        : {number_passing_both}")
    print(f"[done] table written to         : {output_table_path}")


if __name__ == "__main__":
    main()