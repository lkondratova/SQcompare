#!/usr/bin/env python3
import argparse
import pickle
from pathlib import Path
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Import edgeR in R
edgeR = importr("edgeR")

def normalize_expression(expr_dfs):
    """
    Normalize a dictionary of expression dataframes with edgeR TMM.
    Assumes each df has gene_id rows and count columns.
    """
    normalized = {}

    for sample, df in expr_dfs.items():
        if df is None:
            continue

        # First column is assumed to be gene_id
        gene_ids = df.iloc[:, 0]
        counts_only = df.iloc[:, 1:]

        # Use context manager for conversion
        with pandas2ri.converter.context():
            r_counts = pandas2ri.py2rpy(counts_only)
            r_gene_ids = pandas2ri.py2rpy(gene_ids)

        ro.globalenv['counts'] = r_counts
        ro.globalenv['gene_ids'] = r_gene_ids

        # Run edgeR normalization in R
        ro.r('dge <- DGEList(counts=counts)')
        ro.r('dge <- calcNormFactors(dge, method="TMM")')
        norm_r = ro.r('cpm(dge, normalized.lib.sizes=TRUE)')

        # Back-convert to pandas
        with pandas2ri.converter.context():
            norm_df = pandas2ri.rpy2py(norm_r)

        #norm_df.index = gene_ids
        #norm_df.columns = counts_only.columns
        norm_df = pd.DataFrame(norm_df,
                       index=gene_ids,
                       columns=counts_only.columns)

        normalized[sample] = norm_df

    return normalized

def main():
    parser = argparse.ArgumentParser(description="Normalize expression values with edgeR TMM")
    parser.add_argument("--pickle", required=True,
                        help="Pickle file produced by sqanti3_parser.py")
    parser.add_argument("--out", required=True,
                        help="Path to output folder for normalized expression")
    args = parser.parse_args()

    # Load parsed object
    with open(args.pickle, "rb") as f:  # fixed arg name
        parsed = pickle.load(f)

    expr_dfs = {s: parsed["data"][s]["expression"][['universal_id', 'count']] for s in parsed["samples"]}

    # Normalize
    norm_expr = normalize_expression(expr_dfs)

    out_dir = Path(f"{args.out}/normalized_expression")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save per-sample TSVs
    for sample, df in norm_expr.items():
        df.reset_index(inplace=True)
        out_file = out_dir / f"{sample}_normalized_expression.tsv"
        df.to_csv(out_file, sep="\t", index=True)
        parsed["data"][sample]["expression"] = df
        print(f"Saved {out_file}")

    # Save updated pickle with normalized expression
    #parsed["normalized_expression"] = norm_expr
    out_file = f"{args.out}/sqanti3_normalized.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(parsed, f)

    print(f"Normalized expression for {len(norm_expr)} samples")
    print(f"Results saved to {out_file}")

if __name__ == "__main__":
    main()