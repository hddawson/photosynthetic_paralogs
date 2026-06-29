#!/usr/bin/env python3
"""
Find and download PacBio HiFi genome assemblies for flowering plants
(Magnoliopsida = angiosperms, NCBI taxid 3398) from NCBI.

This script shells out to NCBI's official command-line tools:
    - `datasets`    : queries metadata and downloads genome packages
    - `dataformat`  : flattens the metadata JSON into a simple table
Install both (e.g. conda): conda install -c conda-forge ncbi-datasets-cli

WORKFLOW:
    1. Ask NCBI for metadata on every Magnoliopsida assembly.
    2. Keep only rows whose sequencing-tech / assembly-method text looks like HiFi.
    3. Print the matches. By default it stops there (a dry run) so you can
       eyeball the list before committing to a potentially large download.
       Set DRY_RUN = False to actually download.

CAVEAT: the "sequencing tech" field is free text written by whoever submitted
the assembly. It is inconsistent. Some genuine HiFi assemblies do not contain
the word "HiFi" (they may only say "PacBio Sequel II"). So treat the result as
best-effort, and edit HIFI_KEYWORDS below to widen or narrow the net.
"""

import subprocess
import shutil
import csv
import io
import os

# ----------------------------- configuration ------------------------------- #

# NCBI taxon to search. "Magnoliopsida" is the flowering-plant clade (taxid 3398).
target_taxon = "Magnoliopsida"

# Where to put downloaded genome zip files.
output_directory = "hifi_angiosperm_genomes"

# Substrings (lowercased) that we treat as evidence of a HiFi assembly.
# Matched against BOTH the sequencing-tech and assembly-method text.
hifi_keywords = ["hifi", "ccs"]

# If True: only find and print matches, download nothing.
# Flip to False once you are happy with the list.
dry_run = False

# Your NCBI API key. Get one at: https://www.ncbi.nlm.nih.gov/account/
# Doubles the rate limit from 5 to 10 requests/second.
ncbi_api_key = "a99283dd3e3246647ae7116843369336a007"

# --------------------------------------------------------------------------- #


def run_command_capturing_output(command_as_list):
    """Run a shell command and return its stdout as text. Raises on failure."""
    completed = subprocess.run(
        command_as_list, capture_output=True, text=True
    )
    # Surface NCBI's own error message if the call failed, rather than a blank one.
    assert completed.returncode == 0, (
        f"Command failed: {' '.join(command_as_list)}\n{completed.stderr}"
    )
    return completed.stdout


def get_assembly_metadata_table(taxon):
    """
    Return a list of dicts, one per assembly, with the metadata fields we need.

    We pipe `datasets summary ... --as-json-lines` into `dataformat tsv`
    rather than parsing nested JSON ourselves. The requested fields are:
        accession            -> assembly accession (e.g. GCA_xxxxxxxxx.1)
        organism-name        -> species name
        assminfo-release-date     -> assembly release date (used to pick newest)
        assminfo-sequencing-tech  -> free-text sequencing technology
        assminfo-assembly-method  -> free-text assembler / method
    """
    requested_fields = (
        "accession,organism-name,assminfo-release-date,"
        "assminfo-sequencing-tech,assminfo-assembly-method"
    )

    # Step 1: get metadata as one JSON object per line.
    metadata_json_lines = run_command_capturing_output(
        ["datasets", "summary", "genome", "taxon", taxon, "--as-json-lines"]
    )

    # Step 2: flatten to a TSV table via dataformat (reads JSON lines on stdin).
    dataformat_process = subprocess.run(
        ["dataformat", "tsv", "genome", "--fields", requested_fields],
        input=metadata_json_lines,
        capture_output=True,
        text=True,
    )
    assert dataformat_process.returncode == 0, (
        f"dataformat failed:\n{dataformat_process.stderr}"
    )

    # Parse the TSV into dict rows keyed by the human-readable column headers.
    table_rows = list(
        csv.DictReader(io.StringIO(dataformat_process.stdout), delimiter="\t")
    )
    assert len(table_rows) > 0, "NCBI returned no assemblies for this taxon."
    return table_rows


def row_looks_like_hifi(assembly_row):
    """True if any HiFi keyword appears in the seq-tech or assembly-method text."""
    # Column names come from dataformat's headers; join the two relevant ones.
    text_to_search = " ".join(
        assembly_row.get(column_name, "")
        for column_name in ("Assembly Sequencing Tech", "Assembly Method")
    ).lower()
    return any(keyword in text_to_search for keyword in hifi_keywords)


def species_binomial(organism_name):
    """
    Reduce an organism name to its binomial (genus + species epithet) so that
    subspecies and cultivars collapse onto one species.
    e.g. "Zea mays subsp. mays" -> "zea mays".
    To instead keep subspecies separate, return organism_name.strip().lower().
    """
    name_words = organism_name.strip().split()
    return " ".join(name_words[:2]).lower()


def keep_newest_assembly_per_species(assembly_rows):
    """
    Given HiFi assembly rows, return one row per species: the one with the
    latest release date. Release dates are ISO strings (YYYY-MM-DD), which sort
    correctly as plain text, so we compare them directly without date parsing.
    """
    newest_row_by_species = {}
    for current_row in assembly_rows:
        species = species_binomial(current_row["Organism Name"])
        current_date = current_row.get("Assembly Release Date", "")
        existing_row = newest_row_by_species.get(species)
        # Keep current_row if this species is unseen, or if it is more recent.
        if existing_row is None or current_date > existing_row.get(
            "Assembly Release Date", ""
        ):
            newest_row_by_species[species] = current_row
    return list(newest_row_by_species.values())


def download_single_genome(accession, destination_directory):
    """Download one assembly's genome FASTA package as a zip named by accession."""
    output_zip_path = os.path.join(destination_directory, f"{accession}.zip")
    run_command_capturing_output(
        [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",            # FASTA sequence only, no annotation
            "--filename", output_zip_path,
        ]
    )
    return output_zip_path


def main():
    # Fail early with a clear message if the NCBI tools are not installed.
    assert shutil.which("datasets") is not None, \
        "The 'datasets' CLI is not on your PATH. Install ncbi-datasets-cli."
    assert shutil.which("dataformat") is not None, \
        "The 'dataformat' CLI is not on your PATH. Install ncbi-datasets-cli."

    all_assemblies = get_assembly_metadata_table(target_taxon)
    print(f"Retrieved metadata for {len(all_assemblies)} {target_taxon} assemblies.")

    hifi_assemblies = [row for row in all_assemblies if row_looks_like_hifi(row)]
    print(f"Of those, {len(hifi_assemblies)} look like HiFi assemblies.")

    # Collapse to a single newest assembly per species.
    representative_assemblies = keep_newest_assembly_per_species(hifi_assemblies)
    assert len(representative_assemblies) <= len(hifi_assemblies), \
        "Dedup should never increase the number of assemblies."
    print(
        f"After keeping one (newest) per species: "
        f"{len(representative_assemblies)} assemblies.\n"
    )
    for row in sorted(representative_assemblies, key=lambda r: r["Organism Name"]):
        print(
            f"  {row['Assembly Accession']:20s} {row['Organism Name']:35s} "
            f"| {row.get('Assembly Release Date', '')} "
            f"| {row.get('Assembly Sequencing Tech', '')}"
        )

    if dry_run:
        print("\nDRY RUN: nothing downloaded. Set dry_run = False to download.")
        return

    os.makedirs(output_directory, exist_ok=True)
    for row in representative_assemblies:
        accession = row["Assembly Accession"]
        print(f"Downloading {accession} ...")
        download_single_genome(accession, output_directory)
    print(f"\nDone. Genomes saved as zip files in: {output_directory}/")


if __name__ == "__main__":
    main()