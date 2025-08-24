#!/usr/bin/env python3
import argparse
import pickle
from pathlib import Path
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# Activate automatic pandas <-> R dataframe conversion
pandas2ri.activate()

# Import edgeR in R
ro.r('library(edgeR)')

def normalize_expression(expr_dfs):
    """
    Normalize a dictionary of expression dataframes with edgeR TMM.
    Assumes each df has gene_id rows and count/TPM columns.
    """
    normalized = {}

    for sample, df in expr_dfs.items():
        if df is None:
            continue

        # First column is assumed to be gene_id
        gene_ids = df.iloc[:, 0]
        counts_only = df.iloc[:, 1:]

        # R expects a matrix of counts
        r_df = pandas2ri.py2rpy(counts_only)
        ro.globalenv['counts'] = r_df
        ro.globalenv['gene_ids'] = pandas2ri.py2rpy(gene_ids)

        # Run edgeR normalization
        ro.r('dge <- DGEList(counts=counts)')
        ro.r('dge <- calcNormFactors(dge, method="TMM")')
        norm_r = ro.r('cpm(dge, normalized.lib.sizes=TRUE)')

        norm_df = pandas2ri.rpy2py(norm_r)
        norm_df.index = gene_ids
        norm_df.columns = counts_only.columns

        normalized[sample] = norm_df

    return normalized

def main():
    parser = argparse.ArgumentParser(description="Normalize expression values with edgeR TMM")
    parser.add_argument("--input_pickle", required=True,
                        help="Pickle file produced by sqanti3_parser.py")
    parser.add_argument("--out", required=True,
                        help="Path to output folder for normalized expression")
    args = parser.parse_args()

    # Load parsed object
    with open(args.input_pickle, "rb") as f:
        parsed = pickle.load(f)

    expr_dfs = {s: parsed["data"][s]["expression"] for s in parsed["samples"]}

    # Normalize
    norm_expr = normalize_expression(expr_dfs)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save per-sample TSVs
    for sample, df in norm_expr.items():
        out_file = out_dir / f"{sample}_normalized_expression.tsv"
        df.to_csv(out_file, sep="\t", index=True)
        print(f"[INFO] Saved {out_file}")

    # Save updated pickle with normalized expression
    parsed["normalized_expression"] = norm_expr
    out_file = out_dir / "sqanti3_normalized.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(parsed, f)

    print(f"[INFO] Normalized expression for {len(norm_expr)} samples")
    print(f"[INFO] Results saved to {out_file}")

if __name__ == "__main__":
    main()