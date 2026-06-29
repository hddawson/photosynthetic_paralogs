#!/usr/bin/env python3

import subprocess
import shutil
import csv
import io
import os

target_taxon = "Magnoliopsida"
output_directory = "hifi_angiosperm_metadata"

hifi_keywords = ["hifi", "ccs"]

# Do not hard-code API keys in scripts.
# Prefer setting it in your shell:
# export NCBI_API_KEY="your_key_here"
ncbi_api_key = os.environ.get("NCBI_API_KEY", "")


def run_command_capturing_output(command_as_list):
    completed = subprocess.run(
        command_as_list,
        capture_output=True,
        text=True,
        env={**os.environ, "NCBI_API_KEY": ncbi_api_key} if ncbi_api_key else None,
    )
    assert completed.returncode == 0, (
        f"Command failed: {' '.join(command_as_list)}\n{completed.stderr}"
    )
    return completed.stdout


def get_assembly_metadata_table(taxon):
    requested_fields = (
        "accession,organism-name,assminfo-release-date,"
        "assminfo-sequencing-tech,assminfo-assembly-method"
    )

    metadata_json_lines = run_command_capturing_output(
        ["datasets", "summary", "genome", "taxon", taxon, "--as-json-lines"]
    )

    dataformat_process = subprocess.run(
        ["dataformat", "tsv", "genome", "--fields", requested_fields],
        input=metadata_json_lines,
        capture_output=True,
        text=True,
    )
    assert dataformat_process.returncode == 0, (
        f"dataformat failed:\n{dataformat_process.stderr}"
    )

    table_rows = list(
        csv.DictReader(io.StringIO(dataformat_process.stdout), delimiter="\t")
    )
    assert len(table_rows) > 0, "NCBI returned no assemblies for this taxon."
    return table_rows


def row_looks_like_hifi(assembly_row):
    text_to_search = " ".join(
        assembly_row.get(column_name, "")
        for column_name in ("Assembly Sequencing Tech", "Assembly Method")
    ).lower()

    return any(keyword in text_to_search for keyword in hifi_keywords)


def species_binomial(organism_name):
    name_words = organism_name.strip().split()
    return " ".join(name_words[:2]).lower()


def keep_newest_assembly_per_species(assembly_rows):
    newest_row_by_species = {}

    for current_row in assembly_rows:
        species = species_binomial(current_row["Organism Name"])
        current_date = current_row.get("Assembly Release Date", "")
        existing_row = newest_row_by_species.get(species)

        if existing_row is None or current_date > existing_row.get(
            "Assembly Release Date", ""
        ):
            newest_row_by_species[species] = current_row

    return list(newest_row_by_species.values())


def save_csv(rows, output_path):
    assert rows, f"No rows to save for {output_path}"

    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    assert shutil.which("datasets") is not None, \
        "The 'datasets' CLI is not on your PATH. Install ncbi-datasets-cli."
    assert shutil.which("dataformat") is not None, \
        "The 'dataformat' CLI is not on your PATH. Install ncbi-datasets-cli."

    os.makedirs(output_directory, exist_ok=True)

    all_assemblies = get_assembly_metadata_table(target_taxon)
    print(f"Retrieved metadata for {len(all_assemblies)} {target_taxon} assemblies.")

    hifi_assemblies = [row for row in all_assemblies if row_looks_like_hifi(row)]
    print(f"Of those, {len(hifi_assemblies)} look like HiFi assemblies.")

    representative_assemblies = keep_newest_assembly_per_species(hifi_assemblies)
    print(
        f"After keeping one newest assembly per species: "
        f"{len(representative_assemblies)} assemblies."
    )

    all_metadata_path = os.path.join(output_directory, "all_magnoliopsida_assemblies.csv")
    hifi_metadata_path = os.path.join(output_directory, "hifi_representative_assemblies.csv")

    all_metadata_path = os.path.join(
        output_directory,
        "all_magnoliopsida_assemblies.csv"
    )

    all_hifi_path = os.path.join(
        output_directory,
        "all_hifi_like_assemblies.csv"
    )

    representative_hifi_path = os.path.join(
        output_directory,
        "hifi_representative_newest_per_species.csv"
    )

    print(f"All assemblies: {len(all_assemblies)}")
    print(f"All HiFi-like assemblies: {len(hifi_assemblies)}")
    print(f"Newest-per-species HiFi assemblies: {len(representative_assemblies)}")
    print(f"Removed by species deduplication: {len(hifi_assemblies) - len(representative_assemblies)}")

    save_csv(all_assemblies, all_metadata_path)
    save_csv(hifi_assemblies, all_hifi_path)
    save_csv(representative_assemblies, representative_hifi_path)

    print(f"\nSaved all assembly metadata to: {all_metadata_path}")
    print(f"Saved all HiFi-like assemblies to: {all_hifi_path}")
    print(f"Saved newest-per-species HiFi assemblies to: {representative_hifi_path}")
    print("\nNo genomes were downloaded.")

    print(f"\nSaved all assembly metadata to: {all_metadata_path}")
    print(f"Saved HiFi representative metadata to: {hifi_metadata_path}")
    print("\nNo genomes were downloaded.")


if __name__ == "__main__":
    main()