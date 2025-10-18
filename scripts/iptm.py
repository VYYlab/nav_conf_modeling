from tqdm import tqdm
import pandas as pd
from pathlib import Path
from tap import tapify
from dataclasses import dataclass
import json

from tqdm.contrib.concurrent import process_map
tqdm.pandas()

@dataclass
class Arguments:
    source: Path

def process_single(file: Path):
    row = {
        "filename": file.name,
        'group': file.parent.name,
    }

    json_object = json.loads(file.read_text())

    if "iptm" in json_object:
        row["iptm"] = json_object["iptm"]
        return row

def main(args: Arguments):
    source = Path(args.source)
    source = list(source.glob("**/*.json")) if source.is_dir() else [source]

    rows = process_map(process_single, source, chunksize=10)
    rows = [r for r in rows if r]
    df = pd.DataFrame(rows)

    df.to_csv(source[0].parent / "iptm.csv", index=False)

if __name__ == "__main__":
    args = tapify(Arguments)
    main(args)