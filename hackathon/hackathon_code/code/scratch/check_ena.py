#!/usr/bin/env python3
"""
Compare HiFi angiosperm genome assemblies between ENA and NCBI.

BACKGROUND:
ENA (European Nucleotide Archive at EMBL-EBI) and NCBI are both members
of INSDC — they routinely mirror each other's data. However, there can be
a lag, or metadata may differ enough that our keyword filter catches
assemblies in one database that it misses in the other.

This script:
    1. Queries ENA's Portal API for all Magnoliopsida assemblies and
       filters for HiFi keywords in the sequencing_method field.
    2. Queries NCBI via the `datasets` CLI (same logic as the download script).
    3. Compares the two sets by GCA accession and by species name.
    4. Reports:
         - Species found in ENA but not NCBI (potential new downloads)
         - Species found in NCBI but not ENA
         - Accession-level mismatches (same species, different assembly)

REQUIREMENTS:
    - Python standard library only for ENA (uses requests via urllib).
    - NCBI `datasets` and `dataformat` CLIs on PATH for the NCBI side.
      Install: conda install -c conda-forge ncbi-datasets-cli

OUTPUT:
    Prints a comparison report to stdout and writes three TSV files:
        - only_in_ena.tsv
        - only_in_ncbi.tsv
        - in_both_different_accession.tsv
"""

import subprocess
import urllib.request
import urllib.parse
import csv
import io
import json
import os

# ----------------------------- configuration ------------------------------- #

# Magnoliopsida NCBI taxid. Used in the ENA tax_tree() query.
magnoliopsida_taxid = 3398

# Same keywords used in the download script. Keep these in sync.
hifi_keywords = ["hifi", "ccs"]

# Where to write the three output TSV comparison files.
output_directory = "."

# Optional NCBI API key (same as download script). None = no key.
ncbi_api_key = None

# --------------------------------------------------------------------------- #

# ---- ENA Portal API --------------------------------------------------------

# ENA Portal API base URL for assembly searches.
ena_portal_api_base_url = "https://www.ebi.ac.uk/ena/portal/api/search"

# Fields to request from ENA for each assembly record.
# analysis_accession -> ENA analysis accession (ERZ...); GCA not directly available
# scientific_name -> organism name
# last_updated    -> used for recency-based dedup (same logic as NCBI script)
# sequencing_method -> free-text field analogous to NCBI's sequencing_tech
# assembly_software -> assembler used; sometimes mentions HiFi/CCS
ena_fields_to_request = (
    "analysis_accession,scientific_name,last_updated,"
    "sequencing_method,assembly_software"
)


def fetch_ena_assembly_metadata(taxid):
    """
    Query ENA Portal API for all genome assemblies under a given NCBI taxon ID.
    Returns a list of dicts, one per assembly row.

    Uses result=analysis with analysis_type=SEQUENCE_ASSEMBLY, which is the
    correct ENA result type for genome assemblies (result=assembly gives 400).
    The accession field in this context is analysis_accession (e.g. ERZ...).
    GCA_ accessions (INSDC-shared) are assigned later by ENA and may not be
    directly returned here, so we compare by species name rather than accession.

    ENA paginates results; we loop with offset until all records are fetched.
    """
    all_rows = []
    page_size = 10000   # maximum ENA allows per request
    offset = 0

    while True:
        query_parameters = {
            "result": "analysis",
            "query": f'tax_tree({taxid}) AND analysis_type="SEQUENCE_ASSEMBLY"',
            "fields": ena_fields_to_request,
            "format": "tsv",
            "limit": page_size,
            "offset": offset,
        }
        url_with_params = (
            ena_portal_api_base_url + "?" + urllib.parse.urlencode(query_parameters)
        )
        with urllib.request.urlopen(url_with_params) as http_response:
            assert http_response.status == 200, (
                f"ENA Portal API returned status {http_response.status}"
            )
            response_text = http_response.read().decode("utf-8")

        rows_on_this_page = list(
            csv.DictReader(io.StringIO(response_text), delimiter="\t")
        )
        all_rows.extend(rows_on_this_page)
        print(rows_on_this_page)

        # If we got fewer rows than the page size, we have reached the end.
        if len(rows_on_this_page) < page_size:
            break
        offset += page_size

    return all_rows


# ---- NCBI CLI --------------------------------------------------------------


def api_key_flags():
    """Return ['--api-key', key] if one is configured, otherwise empty list."""
    if ncbi_api_key is not None:
        return ["--api-key", ncbi_api_key]
    return []


def run_command_capturing_output(command_as_list):
    """Run a shell command and return its stdout as text. Raises on failure."""
    completed = subprocess.run(command_as_list, capture_output=True, text=True)
    assert completed.returncode == 0, (
        f"Command failed: {' '.join(command_as_list)}\n{completed.stderr}"
    )
    return completed.stdout


def fetch_ncbi_assembly_metadata(taxon_name):
    """
    Query NCBI via the datasets CLI for all assemblies under a taxon name.
    Returns a list of dicts with the same logical fields as the ENA rows.
    """
    requested_fields = (
        "accession,organism-name,assminfo-release-date,"
        "assminfo-sequencing-tech,assminfo-assembly-method"
    )
    metadata_json_lines = run_command_capturing_output(
        ["datasets", "summary", "genome", "taxon", taxon_name, "--as-json-lines"]
        + api_key_flags()
    )
    dataformat_result = subprocess.run(
        ["dataformat", "tsv", "genome", "--fields", requested_fields],
        input=metadata_json_lines,
        capture_output=True,
        text=True,
    )
    assert dataformat_result.returncode == 0, (
        f"dataformat failed:\n{dataformat_result.stderr}"
    )
    rows = list(
        csv.DictReader(io.StringIO(dataformat_result.stdout), delimiter="\t")
    )
    assert len(rows) > 0, "NCBI returned no assemblies."
    return rows


# ---- Shared filtering & normalization --------------------------------------


def row_looks_like_hifi_ena(ena_row):
    """True if HiFi keywords appear in ENA's sequencing_method or assembly_software."""
    text_to_search = " ".join([
        ena_row.get("sequencing_method", ""),
        ena_row.get("assembly_software", ""),
    ]).lower()
    return any(keyword in text_to_search for keyword in hifi_keywords)


def row_looks_like_hifi_ncbi(ncbi_row):
    """True if HiFi keywords appear in NCBI's sequencing tech or assembly method."""
    text_to_search = " ".join([
        ncbi_row.get("Assembly Sequencing Tech", ""),
        ncbi_row.get("Assembly Method", ""),
    ]).lower()
    return any(keyword in text_to_search for keyword in hifi_keywords)


def species_binomial(organism_name):
    """
    Reduce organism name to binomial (first two words), lowercased.
    e.g. "Zea mays subsp. mays" -> "zea mays"
    """
    name_words = organism_name.strip().split()
    return " ".join(name_words[:2]).lower()


def keep_newest_per_species_ena(assembly_rows):
    """
    One assembly per species from ENA: keep the one with the latest last_updated.
    ISO date strings sort correctly as plain text.
    """
    newest_by_species = {}
    for row in assembly_rows:
        species = species_binomial(row.get("scientific_name", ""))
        current_date = row.get("last_updated", "")
        existing = newest_by_species.get(species)
        if existing is None or current_date > existing.get("last_updated", ""):
            newest_by_species[species] = row
    return list(newest_by_species.values())


def keep_newest_per_species_ncbi(assembly_rows):
    """
    One assembly per species from NCBI: keep the one with the latest release date.
    """
    newest_by_species = {}
    for row in assembly_rows:
        species = species_binomial(row.get("Organism Name", ""))
        current_date = row.get("Assembly Release Date", "")
        existing = newest_by_species.get(species)
        if existing is None or current_date > existing.get("Assembly Release Date", ""):
            newest_by_species[species] = row
    return list(newest_by_species.values())


# ---- Comparison ------------------------------------------------------------


def compare_by_species(ena_representatives, ncbi_representatives):
    """
    Compare the two lists by species binomial and GCA accession.

    Returns three lists of dicts:
        only_in_ena       : species with HiFi assembly in ENA, absent from NCBI list
        only_in_ncbi      : species with HiFi assembly in NCBI, absent from ENA list
        different_accession : same species in both, but different GCA accessions
    """
    # Build lookup dicts keyed by species binomial.
    ena_by_species = {
        species_binomial(row["scientific_name"]): row
        for row in ena_representatives
    }
    ncbi_by_species = {
        species_binomial(row["Organism Name"]): row
        for row in ncbi_representatives
    }

    all_species = set(ena_by_species) | set(ncbi_by_species)

    only_in_ena = []
    only_in_ncbi = []
    different_accession = []

    for species in sorted(all_species):
        in_ena = species in ena_by_species
        in_ncbi = species in ncbi_by_species

        if in_ena and not in_ncbi:
            only_in_ena.append(ena_by_species[species])
        elif in_ncbi and not in_ena:
            only_in_ncbi.append(ncbi_by_species[species])
        else:
            # Species appears in both — check if the accessions match.
            ena_accession = ena_by_species[species].get("analysis_accession", "")
            ncbi_accession = ncbi_by_species[species].get("Assembly Accession", "")
            # GCA accessions from both databases should be identical (INSDC sharing).
            # If they differ it may mean one database has a newer assembly version.
            if ena_accession != ncbi_accession:
                different_accession.append({
                    "species": species,
                    "ena_accession": ena_accession,
                    "ena_last_updated": ena_by_species[species].get("last_updated", ""),
                    "ncbi_accession": ncbi_accession,
                    "ncbi_release_date": ncbi_by_species[species].get(
                        "Assembly Release Date", ""
                    ),
                })

    return only_in_ena, only_in_ncbi, different_accession


# ---- TSV writing -----------------------------------------------------------


def write_tsv(rows_as_dicts, output_path):
    """Write a list of dicts to a TSV file."""
    if not rows_as_dicts:
        print(f"  (nothing to write for {output_path})")
        return
    with open(output_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(
            tsv_file, fieldnames=rows_as_dicts[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows_as_dicts)
    print(f"  Written: {output_path}")


# ---- Main ------------------------------------------------------------------


def main():
    import shutil as _shutil
    assert _shutil.which("datasets") is not None, \
        "The 'datasets' CLI is not on your PATH. Install ncbi-datasets-cli."
    assert _shutil.which("dataformat") is not None, \
        "The 'dataformat' CLI is not on your PATH. Install ncbi-datasets-cli."

    # --- ENA ---
    print(f"Querying ENA for Magnoliopsida assemblies (taxid {magnoliopsida_taxid})...")
    all_ena_rows = fetch_ena_assembly_metadata(magnoliopsida_taxid)
    print(f"  ENA returned {len(all_ena_rows)} total assemblies.")

    hifi_ena_rows = [row for row in all_ena_rows if row_looks_like_hifi_ena(row)]
    print(f"  Of those, {len(hifi_ena_rows)} look like HiFi.")

    ena_representatives = keep_newest_per_species_ena(hifi_ena_rows)
    print(f"  After dedup (one per species): {len(ena_representatives)}.")

    # --- NCBI ---
    print("\nQuerying NCBI for Magnoliopsida assemblies...")
    all_ncbi_rows = fetch_ncbi_assembly_metadata("Magnoliopsida")
    print(f"  NCBI returned {len(all_ncbi_rows)} total assemblies.")

    hifi_ncbi_rows = [row for row in all_ncbi_rows if row_looks_like_hifi_ncbi(row)]
    print(f"  Of those, {len(hifi_ncbi_rows)} look like HiFi.")

    ncbi_representatives = keep_newest_per_species_ncbi(hifi_ncbi_rows)
    print(f"  After dedup (one per species): {len(ncbi_representatives)}.")

    # --- Compare ---
    print("\nComparing by species binomial...")
    only_in_ena, only_in_ncbi, different_accession = compare_by_species(
        ena_representatives, ncbi_representatives
    )

    print(f"\n{'='*60}")
    print(f"  Species in ENA only  (not in NCBI HiFi list): {len(only_in_ena)}")
    print(f"  Species in NCBI only (not in ENA HiFi list):  {len(only_in_ncbi)}")
    print(f"  Same species, different accession:            {len(different_accession)}")
    print(f"{'='*60}\n")

    if only_in_ena:
        print("--- ENA only (potential missing from your NCBI download) ---")
        for row in only_in_ena:
            print(
                f"  {row['analysis_accession']:20s} {row['scientific_name']:35s} "
                f"| {row.get('sequencing_method', '')}"
            )

    if only_in_ncbi:
        print("\n--- NCBI only (not yet in ENA, or missing ENA sequencing_method) ---")
        for row in only_in_ncbi:
            print(
                f"  {row['Assembly Accession']:20s} {row['Organism Name']:35s} "
                f"| {row.get('Assembly Sequencing Tech', '')}"
            )

    if different_accession:
        print("\n--- Same species, different accession (one DB may be newer) ---")
        for row in different_accession:
            print(
                f"  {row['species']:35s}  ENA: {row['ena_accession']} ({row['ena_last_updated']}) "
                f" NCBI: {row['ncbi_accession']} ({row['ncbi_release_date']})"
            )

    # Write TSV outputs.
    print("\nWriting TSV files...")
    write_tsv(only_in_ena,         os.path.join(output_directory, "only_in_ena.tsv"))
    write_tsv(only_in_ncbi,        os.path.join(output_directory, "only_in_ncbi.tsv"))
    write_tsv(different_accession, os.path.join(output_directory, "in_both_different_accession.tsv"))


if __name__ == "__main__":
    main()
