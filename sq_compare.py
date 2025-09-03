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
    
    # Make output folders
    os.makedirs(args.out, exist_ok=True)
    
    # Step 1: Parse inputs
    os.makedirs(parsed_dir, exist_ok=True)
    run_script("scripts/parse_sq_inputs.py", ["--input_files", args.input_files, "--out", args.out])
    
    # Step 2: Collapse ISM (optional)
    if args.collapseISM:
        run_script("collapse_ism.py", ["--pickle", f'{args.out}/sqanti3_samples.pkl', "--out", args.out])
    
    # Step 3: Assign universal IDs
    if args.collapseISM:
        pickle_df = f'{args.out}/sqanti3_samples_ISMcollapsed.pkl' # use collapsed pickle if ISM collapsing was done
    else:
        pickle_df = f'{args.out}/sqanti3_samples.pkl'
 
    run_script("scripts/universal_id.py", ["--pickle", pickle_df, "--out", args.out])
    
    # Step 4: TMM normalization (optional)
    if args.expression_provided:
        os.makedirs(norm_dir, exist_ok=True)
        pickle_df = f'{args.out}/sqanti3_standardizeds.pkl'
        run_script("tmm_norm.py", ["--input_files", pickle_df, "--out", args.out])
    





    
    
    # Step 5: Create matrix and isoform info
    gen_dir = os.path.join(args.outdir, "generalized")
    os.makedirs(gen_dir, exist_ok=True)
    run_script("generalize_isoforms.py", ["--input_files", working_dir, "--out", gen_dir])
    
    # Step 6: Generate plots and summary
    plots_dir = os.path.join(args.outdir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    run_script("sq_compare_summary.py", ["--input_files", gen_dir, "--out", plots_dir])
    
    # Step 7: Export final report
    export_dir = os.path.join(args.outdir, "final_report")
    os.makedirs(export_dir, exist_ok=True)
    run_script("export_script.py", ["--input", plots_dir, "--out", export_dir])
    
    print(f"Pipeline finished! Results are in: {args.outdir}")

if __name__ == "__main__":
    main()

#python main.py --input_files inputs.tsv --out results --collapseISM