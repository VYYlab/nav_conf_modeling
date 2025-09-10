# python3 subset/subset_plddt.py --source "/Users/kushnarang/Downloads/nav19models/_hnav1.9_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb" --metrics subset/sections.txt

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
import statistics

class Arguments(Tap):
    source: str
    metrics: str

def parse_selector(metric):
    rrange = metric["range"].replace(" ", "")
    sections = rrange.split(",")
    all_residues = []
    for section in sections:
        start, end = section.split("-")
        start = int(start[1:])
        end = int(end[1:])
        all_residues.extend(list(range(start, end+1)))
    
    metric["residues"] = all_residues
    return metric

def process_single(processing_arguments):
    file, selectors = processing_arguments

    row = {
        "filename": file.name,
        'group': file.parent.name,
    }

    pose = pose_from_pdb(str(file))
    info = pose.pdb_info()

    for metric in selectors:
        tag = metric["tag"]
        plddts = []
        for residue in metric["residues"]:
            plddts.append(info.bfactor(residue, 1))
        row[tag] = statistics.fmean(plddts)
    
    return row

def main():
    args = Arguments().parse_args()
    source = Path(args.source)
    source = list(source.glob("**/*.pdb")) if source.is_dir() else [source]

    metrics = pd.read_csv(args.metrics, delim_whitespace=True).to_dict(orient='records')
    selectors = [parse_selector(metric) for metric in metrics]

    processing_arguments = [(file, selectors) for file in source]
    rows = process_map(process_single, processing_arguments, chunksize=10)

    df = pd.DataFrame(rows)
    df.to_csv(source[0].parent / "subset_plddts_outer.csv", index=False)

if __name__ == "__main__":
    main()

