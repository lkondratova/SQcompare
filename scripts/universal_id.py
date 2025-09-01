#!/usr/bin/env python3
import argparse
import pickle
from pathlib import Path
import pandas as pd

"""
python standardize_isoform_ids_gtf.py \
    --pickle results/sqanti3_normalized.pkl \
    --out output_folder

"""

def extract_junction_chain_from_gtf(gtf_file):
    """
    Extract junction chains from a GTF file.
    Returns dict: transcript_id -> list of junction coordinates [chr, start1, end1, start2, end2,...]
    """
    chains = {}
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None,
                      names=["chr", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attribute"])
    
    exons = gtf[gtf["feature"] == "exon"].copy()
    
    # Extract transcript_id
    def get_transcript_id(attr):
        for field in attr.split(";"):
            if 'transcript_id' in field:
                return field.split('"')[1]
        return None
    
    exons["transcript_id"] = exons["attribute"].apply(get_transcript_id)
    
    # Group by transcript_id
    for tid, df in exons.groupby("transcript_id"):
        df_sorted = df.sort_values("start")
        # Junction chain: [chr, start1, end1, start2, end2, ...]
        chain = [df_sorted.iloc[0]["chr"]]
        for _, row in df_sorted.iterrows():
            chain.extend([row["start"], row["end"]])
        chains[tid] = chain

    return chains

def standardize_isoforms_cross_sample(pickle, out_dir=None):
    """
    Standardize isoform IDs across all samples using junction chains from GTF files.
    """
    samples = pickle["data"]["samples"]
    # extract junction chains per sample
    sample_chains = {}
    for sample in samples:
        sample_chains[sample] = extract_junction_chain_from_gtf(pickle["data"][sample]["gtf"])
        print(f"Extracted {len(sample_chains[sample])} junction chains for sample {sample}")
    # collect all junction chains into a set of unique chains
    all_chains_set = set()
    for chains in sample_chains.values():
        for chain in chains.values():
            all_chains_set.add(tuple(chain))  # convert list -> tuple for hashability

    # assign universal IDs
    unique_chains = list(all_chains_set)
    ids = {tuple(chain): f"isoform{i+1}" for i, chain in enumerate(unique_chains)}

    # create sample-specific mapping and add to DataFrames
    for sample in samples:
        chains = sample_chains[sample]
        iso_map = {tid: ids[tuple(chain)] for tid, chain in chains.items()}

        # Update junctions
        junc_df = pickle["data"][sample]["junctions"].copy()
        junc_df["universal_id"] = junc_df["isoform"].map(iso_map)
        pickle["data"][sample]["junctions"] = junc_df

        # Update classification
        class_df = pickle["data"][sample]["classification"].copy()
        class_df["universal_id"] = class_df["isoform"].map(iso_map)
        pickle["data"][sample]["classification"] = class_df

        # Update expression if available
        expr_df = pickle["data"][sample]["expression"]
        if expr_df is not None:
            expr_df = expr_df.copy()
            expr_df.insert(0, "universal_id", expr_df.iloc[:, 0].map(iso_map))
            pickle["data"][sample]["expression"] = expr_df

        # save TSVs
        junc_df.to_csv(out_dir / f"{sample}_junctions_std.tsv", sep="\t", index=False)
        class_df.to_csv(out_dir / f"{sample}_classification_std.tsv", sep="\t", index=False)
        #if expr_df is not None:
        #    expr_df.to_csv(out_dir / f"{sample}_expression_std.tsv", sep="\t", index=False)

    return pickle

def main():
    parser = argparse.ArgumentParser(description="Standardize isoform IDs across all samples using GTF")
    parser.add_argument("--pickle", required=True,
                        help="Pickle file with parsed SQANTI3 outputs")
    parser.add_argument("--out", required=True,
                        help="Output folder for updated TSVs and pickle")

    args = parser.parse_args()

    out_dir = Path(args.out)

    # Load parsed data
    with open(args.pickle, "rb") as f:
        parsed_obj = pickle.load(f)
        # Standardize isoforms
        updated_obj = standardize_isoforms_cross_sample(args.pickle, args.out)

    # Save updated pickle
    out_pickle = f"{out_dir}/sqanti3_standardized.pkl"
    with open(args.pickle, "wb") as f:
        pickle.dump(updated_obj, f)

    print(f"[INFO] Standardized isoform IDs for {len(parsed_obj['samples'])} samples")
    print(f"[INFO] Updated pickle saved to {out_pickle}")
    #print(f"[INFO] TSVs saved in {out_dir}")

if __name__ == "__main__":
    main()