#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Create isoform metadata and presence/expression matrices")
    parser.add_argument(
        "--pickle", required=True,
        help="Paths to input pickle with classification, junctions, and expression dataframes"
    )
    parser.add_argument(
        "--out", required=True,
        help="Path to the output folder"
    )
    args = parser.parse_args()

    os.makedirs(f"{args.out}/summarized", exist_ok=True)

    all_classifications = []
    all_expr = {}
    with open(args.pickle, "rb") as f:
        pickle_df = pickle.load(f)
        samples = pickle_df["samples"]
        for sample in samples:
            class_df = pickle_df["data"][sample]["classification"].copy()
            class_df["sample"] = sample
            all_classifications.append(class_df)
            expr_df = pickle_df["data"][sample]["expression"]
            if expr_df is not None:
                expr_df = expr_df[["universal_id", "count"]].copy()
                all_expr[sample] = expr_df.set_index("universal_id")["count"]

        combined = pd.concat(all_classifications, ignore_index=True)
        combined["junction_chain"] = combined["junction_chain"].apply(tuple)

        # --- Isoform info (metadata) ---
        isoform_info = (
            combined
            .groupby("junction_chain")
            .agg({
                "universal_id": "first",   # consistent across samples
                "structural_category": "first",
                "associated_gene": "first",
                "associated_transcript": "first",
                "exons": "first",          # consistent
                "length": "mean"           # averaged
            })
            .reset_index()
            .rename(columns={
                "junction_chain": "unique_jc",
                "structural_category": "category",
                "exons": "exons_n",
                "length": "average_length"
            })
        )

        # --- Isoform matrix (samples x isoforms) ---
        isoform_ids = isoform_info["universal_id"].tolist()
        matrix = pd.DataFrame(index=isoform_ids)

        for sample in combined["sample"].unique():
            if sample in all_expr:
                # expression available
                expr_series = all_expr[sample]
                if expr_series.index.has_duplicates:
                    expr_series = expr_series.groupby(expr_series.index).sum()
                col = expr_series.reindex(isoform_ids).fillna(0)
                matrix[sample] = col
            else:
                # binary presence/absence
                present_ids = set(combined.loc[combined["sample"] == sample, "universal_id"])
                matrix[sample] = [1 if uid in present_ids else 0 for uid in isoform_ids]

        matrix = matrix.reset_index().rename(columns={"index": "universal_id"})

        # --- Save ---
        isoform_info.to_csv(os.path.join(f"{args.out}/summarized", "isoform_info.tsv"), sep="\t", index=False)
        matrix.to_csv(os.path.join(f"{args.out}/summarized", "isoform_matrix.tsv"), sep="\t", index=False)
        with open(os.path.join(args.out, "combined.pkl"), "wb") as f:
            pickle.dump({"isoform_info": isoform_info, "isoform_matrix": matrix}, f)

if __name__ == "__main__":
    main()