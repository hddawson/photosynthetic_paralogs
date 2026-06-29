#!/usr/bin/env python3
import gzip
import re
from pathlib import Path

IDS = Path("data/housekeeping_candidates_Zm.txt")
OUT = Path("data/housekeeping_candidates_Zm.peptides.faa")
FA = Path("data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.gz")

if not FA.exists():
    raise SystemExit(f"Missing {FA}. Download it first with curl.")

wanted = {x.strip().split()[0] for x in IDS.read_text().splitlines() if x.strip()}

hits = {}
current_header = None
current_gene = None
seq = []

def save_record():
    if current_gene in wanted and current_gene not in hits:
        hits[current_gene] = (current_header, "".join(seq))

with gzip.open(FA, "rt") as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            save_record()
            current_header = line[1:]
            seq = []
            m = re.search(r"(Zm00001eb\d+)", current_header)
            current_gene = m.group(1) if m else None
        else:
            seq.append(line)

save_record()

with OUT.open("w") as out:
    for gene in sorted(wanted):
        if gene in hits:
            header, peptide = hits[gene]
            out.write(f">{gene} {header}\n")
            for i in range(0, len(peptide), 60):
                out.write(peptide[i:i+60] + "\n")

missing = sorted(wanted - hits.keys())

print(f"Wrote {len(hits)} peptide sequences to {OUT}")
print(f"Missing: {len(missing)}")

if missing:
    miss_out = Path(str(OUT) + ".missing.txt")
    miss_out.write_text("\n".join(missing) + "\n")
    print(f"Missing IDs saved to {miss_out}")