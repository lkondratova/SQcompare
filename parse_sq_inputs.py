#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import hashlib
import pickle
from pathlib import Path

def get_sample_name(class_file):
    """Extract sample name from classification file (strip extension and suffix)."""
    base = Path(class_file).stem
    if base.endswith("_classification"):
        return base.replace("_classification", "")
    return base

def hash_files(*files):
    """Generate MD5 hash from file paths (concatenated)."""
    m = hashlib.md5()
    for f in files:
        if f is not None:
            m.update(f.encode("utf-8"))
    return m.hexdigest()

#############
def main():
    parser = argparse.ArgumentParser(
        description="Parse multiple SQANTI3 outputs from a TSV mapping file."
    )
    parser.add_argument(
        "--collapse_ISM", type=lambda x: str(x).lower() in ["true", "1", "yes"],
        default=False, help="Collapse ISM isoforms into one category (default: False)"
    )
    parser.add_argument(
        "--input_files", required=True,
        help="TSV file with columns: classification.txt, junctions.txt, [expression.tsv]"
    )
    parser.add_argument(
        "out", required=True
        help="Path to the output folder"
    )
    args = parser.parse_args()

    # Load mapping file
    input_df = pd.read_csv(args.input_files, sep="\t", header=None, comment="#")

    if input_df.shape[1] < 2:
        raise ValueError("TSV must have at least 2 columns: classification.txt, junctions.txt")

    input_df.columns = ["classification", "junctions"] + (
        ["expression"] if input_df.shape[1] >= 3 else []
    )

    samples_info = {}
    for _, row in input_df.iterrows():
        class_file = row["classification"]
        junc_file = row["junctions"]
        expr_file = row["expression"] if "expression" in row else None

        # File existence check
        if not os.path.exists(class_file):
            raise FileNotFoundError(f"Classification file not found: {class_file}")
        if not os.path.exists(junc_file):
            raise FileNotFoundError(f"Junction file not found: {junc_file}")
        if expr_file and not os.path.exists(expr_file):
            raise FileNotFoundError(f"Expression file not found: {expr_file}")

        # Extract sample name
        sample_name = get_sample_name(class_file)

        # Create hash
        file_hash = hash_files(class_file, junc_file, expr_file)

        # Store paths (lazy loading – don’t read files yet, unless needed)
        samples_info[sample_name] = {
            "hash": file_hash,
            "classification_path": class_file,
            "junctions_path": junc_file,
            "expression_path": expr_file
        }

    # Summary
    print(f"Loaded {len(samples_info)} samples")
    print("Sample names:")
    for name in samples_info.keys():
        print(f" - {name}")
    print(f"\n--collapse_ISM = {args.collapse_ISM}")

    # Save intermediate object
    with open("sqanti3_samples.pkl", "wb") as f:
        pickle.dump(samples_info, f)

if __name__ == "__main__":
    main()