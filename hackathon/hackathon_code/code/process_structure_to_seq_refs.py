#!/usr/bin/env python3
"""
Fetch all UniProt protein sequences for a PDB entry and write a combined FASTA.

Usage:
    python fetch_pdb_reference_proteins.py 9LK4

Output:
    data/reference_protein_sequences/<PDB_ID>_reference_proteins.faa
"""

import urllib.request
import json
import os
import time
import sys

# Output directory for reference protein FASTAs.
output_directory = "data/reference_protein_sequences"

rcsb_entry_api_url_template          = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
rcsb_polymer_entity_api_url_template = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
uniprot_fasta_url_template           = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def fetch_uniprot_ids_from_pdb(pdb_id):
    entry_url = rcsb_entry_api_url_template.format(pdb_id=pdb_id.upper())
    with urllib.request.urlopen(entry_url) as http_response:
        entry_data = json.loads(http_response.read().decode("utf-8"))

    polymer_entity_ids = (
        entry_data
        .get("rcsb_entry_container_identifiers", {})
        .get("polymer_entity_ids", [])
    )
    assert len(polymer_entity_ids) > 0, \
        f"No polymer entity IDs found for PDB entry {pdb_id}"

    seen_uniprot_ids    = set()
    ordered_uniprot_ids = []

    for entity_id in polymer_entity_ids:
        entity_url = rcsb_polymer_entity_api_url_template.format(
            pdb_id=pdb_id.upper(), entity_id=entity_id
        )
        with urllib.request.urlopen(entity_url) as http_response:
            entity_data = json.loads(http_response.read().decode("utf-8"))

        uniprot_ids_for_entity = (
            entity_data
            .get("rcsb_polymer_entity_container_identifiers", {})
            .get("uniprot_ids") or []
        )
        for uniprot_id in uniprot_ids_for_entity:
            if uniprot_id not in seen_uniprot_ids:
                seen_uniprot_ids.add(uniprot_id)
                ordered_uniprot_ids.append(uniprot_id)

        time.sleep(0.1)

    return ordered_uniprot_ids


def fetch_uniprot_fasta(uniprot_accession):
    url = uniprot_fasta_url_template.format(accession=uniprot_accession)
    with urllib.request.urlopen(url) as http_response:
        fasta_text = http_response.read().decode("utf-8")
    assert fasta_text.strip().startswith(">"), \
        f"Unexpected response from UniProt for {uniprot_accession}: {fasta_text[:100]}"
    return fasta_text


def main():
    assert len(sys.argv) == 2, "Usage: python fetch_pdb_reference_proteins.py <PDB_ID>"
    pdb_accession     = sys.argv[1].upper()
    output_fasta_path = os.path.join(
        output_directory, f"{pdb_accession}_reference_proteins.faa"
    )
    os.makedirs(output_directory, exist_ok=True)

    print(f"Querying RCSB for UniProt IDs in PDB entry {pdb_accession}...")
    uniprot_ids = fetch_uniprot_ids_from_pdb(pdb_accession)
    assert len(uniprot_ids) > 0, \
        f"No UniProt IDs found for {pdb_accession}."

    print(f"Found {len(uniprot_ids)} unique UniProt accessions:")
    for uniprot_id in uniprot_ids:
        print(f"  {uniprot_id}")

    print("\nDownloading sequences from UniProt...")
    all_fasta_text = []
    for uniprot_id in uniprot_ids:
        print(f"  Fetching {uniprot_id}...")
        all_fasta_text.append(fetch_uniprot_fasta(uniprot_id).strip())
        time.sleep(0.3)

    with open(output_fasta_path, "w") as output_file:
        output_file.write("\n".join(all_fasta_text) + "\n")

    print(f"\nWrote {len(uniprot_ids)} sequences to {output_fasta_path}")


if __name__ == "__main__":
    main()