# Usage:
# runscript rmsd.py full/path/to/models/directory REF1 REF2

# Example:
# runscript rmsd.py ~/Documents/nav1.1-models 7DTD

from collections import defaultdict
from pathlib import Path
from chimerax.core.commands import run
import csv
import sys

directory, *tags = sys.argv[1:]

ACTIVE = directory
ACTIVE_REFERENCES = [{ "source": t, "tag": t } for t in tags]

SOURCE_DIR = ACTIVE
SAVE_NAME = ACTIVE + "/rmsd.csv"

source = Path(SOURCE_DIR)
source = list(source.glob("**/*.pdb")) if source.is_dir() else [source]

data = defaultdict(dict)

for structure in ACTIVE_REFERENCES:
    reference_source = structure["source"]
    tag = structure["tag"]
    run(session, f"open {reference_source}")
    for i, file in enumerate(source):
        if i > 0: run(session, f"del #2")
        run(session, f"open {file}")
        try:
            run(session, "del /B")
        except:
            pass
        x = run(session, f"matchmaker #1 to #2")
        data[str(file)][tag] =  x[0]["final RMSD"]
        data[str(file)][tag + "_full"] =  x[0]["full RMSD"]
        print("-" * 25 + " " + str(i) + " " + "-" * 25)
    run(session, f"del all")

rows = []

for key, value in data.items():
    rows.append({
        "source": key,
        **value
    })

with open(SAVE_NAME, 'w') as csvfile:
    fieldnames = ['source'] + [r["tag"] for r in ACTIVE_REFERENCES] + [r["tag"] + "_full" for r in ACTIVE_REFERENCES]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    print(f"Results saved to {SAVE_NAME}")
