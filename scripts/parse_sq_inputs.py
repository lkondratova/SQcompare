#!/usr/bin/env python3
import argparse
import pandas as pd
import pickle
from pathlib import Path

def parse_sqanti3_inputs(tsv_file):
    """
    Parse a TSV with paths to SQANTI3 outputs.
    
    Parameters
    ----------
    tsv_file : str or Path
        Path to TSV file with 2-3 columns:
        1) classification.txt
        2) junctions.txt
        3) corrected.gtf
        4) (optional) expression levels
    Returns
    -------
    dict
        Dictionary keyed by sample name. Each value is another dict with:
        - classification : DataFrame
        - junctions : DataFrame
        - gtf : corrected.gtf
        - expression : DataFrame or None
    """
    samples_info = {}

    df_inputs = pd.read_csv(tsv_file, sep="\t", header=None)
    n_samples = len(df_inputs)

    for idx, row in df_inputs.iterrows():
        class_path = Path(row[0])
        sj_path = Path(row[1])
        gtf_path = Path(row[2])
        expr_path = Path(row[3]) if len(row) > 3 else None

        # sample name extraction (strip _classification.txt)
        sample_name = class_path.stem.replace("_classification", "")

        samples_info[sample_name] = {
            "classification": pd.read_csv(class_path, sep="\t"),
            "junctions": pd.read_csv(sj_path, sep="\t"),
            "gtf":pd.read_csv(gtf_path, sep="\t", header=None),
            "expression": pd.read_csv(expr_path, sep="\t", header=None) if expr_path else None
        }

    return {
        "n_samples": n_samples,
        "samples": list(samples_info.keys()),
        "data": samples_info
    }

def main():
    parser = argparse.ArgumentParser(description="Parse SQANTI3 outputs from TSV file")
    parser.add_argument("--input_files", required=True,
                        help="TSV file with SQANTI3 output paths")
    parser.add_argument(
        "--out", required=True,
        help="Path to the output folder"
    )
    args = parser.parse_args()

    result = parse_sqanti3_inputs(args.input_files)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "sqanti3_samples.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(result, f)

    print(f"[INFO] Parsed {result['n_samples']} samples: {', '.join(result['samples'])}")
    print(f"[INFO] Results saved in {out_file}")

if __name__ == "__main__":
    main()