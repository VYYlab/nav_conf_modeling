import fnmatch
import pandas as pd
from pathlib import Path
from tap import tapify
from dataclasses import dataclass

@dataclass
class Arguments:
    source: Path

def main(args: Arguments):
    log_file = Path(args.source)
    log_rows = log_file.read_text().splitlines()

    recycle_rows = fnmatch.filter(log_rows, "*alphafold2_multimer_v3_model_*recycle=*")

    def process_row(row_text: str):
        splits = row_text.split(" ")
        model = splits[2]
        iptm = [s for s in splits if "ipTM" in s][0].replace("ipTM=", "")
        recycle = [s for s in splits if "recycle" in s][0].replace("recycle=", "")
        return { "file_id": f"{model}.r{recycle}", "iptm": float(iptm) }

    print(len(recycle_rows))

    rows = list(map(process_row, recycle_rows))
    df = pd.DataFrame(rows)

    df.to_csv(log_file.parent / "iptm_by_log.csv", index=False)

if __name__ == "__main__":
    args = tapify(Arguments)
    main(args)