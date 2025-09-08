#!/usr/bin/env python

import os
import argparse
import subprocess
import sys

def run_script(script, args):
    """Helper function to run a Python script via subprocess."""
    cmd = [sys.executable, script] + args
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, check=True)
    return result

def main():
    parser = argparse.ArgumentParser(description="Run the isoform analysis pipeline without Nextflow.")
    
    parser.add_argument("--input_files", required=True,
                        help="TSV file with paths to SQANTI3 outputs (classification, junctions, GTF, optional expression).")
    parser.add_argument("--out", required=True,
                        help="Output folder for results.")
    parser.add_argument("--collapseISM", action="store_true",
                        help="Whether to collapse ISM isoforms.")
    
    args = parser.parse_args()
    
    #Make output folders
    os.makedirs(args.out, exist_ok=True)
    
    #1: Parse inputs
    run_script("scripts/parse_sq_inputs.py", ["--input_files", args.input_files, "--out", args.out])
    
    #2: Collapse ISM (optional)
    if args.collapseISM:
        run_script("scripts/collapse_ism.py", ["--pickle", f'{args.out}/sqanti3_samples.pkl', "--out", args.out])
    
    #3: Assign universal IDs
    if args.collapseISM:
        pickle_df = f'{args.out}/sqanti3_samples_ISMcollapsed.pkl' # use collapsed pickle if ISM collapsing was done
    else:
        pickle_df = f'{args.out}/sqanti3_samples.pkl'
 
    run_script("scripts/universal_id.py", ["--pickle", pickle_df, "--out", args.out])

    #4 TMM normalization of expression values if provided
    with open(args.input_files) as f:
        first_line = f.readline()
        num_cols = len(first_line.strip().split('\t'))
    if num_cols == 4:
        run_script("scripts/tmm_norm.py", ["--pickle", f"{args.out}/sqanti3_standardized.pkl", "--out", args.out])
        #5 create matrix and isoform info
        run_script("scripts/generalize_isoforms.py", ["--pickle", f"{args.out}/sqanti3_normalized.pkl", "--out", args.out])
    else:
        run_script("scripts/sq_standardized.py", ["--pickle", f"{args.out}/sqanti3_standardized.pkl", "--out", args.out])

    #5: Visualize comparisons
    run_script("scripts/sq_compare_summary.py", ["--out", args.out])

    # Delete all pickle files in the output directory
    for file in os.listdir(args.out):
        if file.endswith(".pkl"):
            os.remove(os.path.join(args.out, file))
if __name__ == "__main__":
    main()