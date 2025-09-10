import pyrosetta
from pyrosetta import pose_from_pdb, rosetta
pyrosetta.init("-mute all")

import numpy as np
from tqdm import tqdm

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from tap import Tap

from tqdm.contrib.concurrent import process_map
tqdm.pandas()

import matplotlib.pyplot as plt

class Arguments(Tap):
    source: str
    metrics: str

def process_single(processing_arguments):
    # score = 
    file, metrics = processing_arguments

    row = {
        "filename": file.name,
        'group': file.parent.name,
    }

    pose = pose_from_pdb(str(file))

    for metric in metrics:
        residue_from = np.array(pose.residue(metric["res1"]).xyz(metric["atom1"]))
        residue_to = np.array(pose.residue(metric["res2"]).xyz(metric["atom2"]))
        distance = np.linalg.norm(residue_to-residue_from)
        row[metric["tag"]] = distance
    
    return row

def main():
    args = Arguments().parse_args()
    source = Path(args.source)
    source = list(source.glob("**/*.pdb")) if source.is_dir() else [source]

    metrics = pd.read_csv(args.metrics, delim_whitespace=True).to_dict(orient='records')

    processing_arguments = [(file, metrics) for file in source]
    rows = process_map(process_single, processing_arguments, chunksize=10)

    df = pd.DataFrame(rows)

    df.to_csv(source[0].parent / "distances.csv", index=False)

    # plt.scatter(df["GC"], df["WIDTH"])

    # plt.title(source[0].parent.name)

    # plt.xlabel("GC")
    # plt.ylabel("WIDTH")

    # plt.show()
    # plt.savefig(source[0].parent / "metrics.png")

if __name__ == "__main__":
    main()

