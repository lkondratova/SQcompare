# SQcompare

This repository contains a Python-based workflow to analyze and compare unique junction chains (UJC) from multiple SQANTI3 outputs. The pipeline handles multiple samples, collapses incomplete splice match (ISM) isoforms (optional), assigns universal isoform IDs across samples, normalizes expression, and generates summary plots and reports.

---

## Features

- Parse SQANTI3 output files: classification, junctions, GTF, expression (optional).
- Collapse ISM isoforms (optional).
- Assign **universal isoform IDs** across samples based on junction chains.
- Normalize expression values using TMM (edgeR-like) normalization if expression files are provided.
- Generate combined isoform matrices and summary files.
- Produce plots and tables:
  - UJCs counts per category
  - Length distributions
  - Heatmap of the top 1000 variable UJCs expression
  - UpSet plots for UJCs sharing
  - Monoexon vs multiexon counts

---

## Installation

This pipeline is designed to run in a **Conda environment**. Use the provided `sq_compare_environment.yml`:

```bash
conda env create -f sq_compare_environment.yml
conda activate sq_compare
```

---
## Input Files

SQcompare accepts a tab-separated input file (no header) with the following columns:

1. Path to a SQANTI3 *classification.txt

2. Path to a SQANTI3 *junctions.txt

3. Path to *corrected.gtf

4. Optionally, an expression file (absolute values) can be provided, where the first column is an isoform name and the second column is an absolute expression value (no header).

An example bash script for generating an input file can be found in `helper_scripts/create_input_file.sh`

---

## Workflow Overview

1. Parse inputs (`parse_sq_inputs.py`). Parse SQANTI3 files and organize them into dataframes.
2. Collapse ISM isoforms, optional (`collapse_ism.py`). Collapse incomplete splice match isoforms if requested. Collapses all FSMs and ISMs of the same transcripts to one of the closest to the reference match; the expression values, if provided, are collapsed accordingly. 
3. Assign universal IDs (`universal_id.py`). Assign universal isoform IDs across all samples based on junction chains, e.g., isoform1, isoform2.
4. Normalize expression if expression values are provided (`tmm_norm.py`). Normalize expression values using TMM edgeR-like normalization.
5. Generate matrices and combined isoform info (`generalize_isoforms.py`). Create combined isoform matrices and information files for all samples.
   - isoform_info.tsv: unique_jc, universal_id, category, associated_gene, associated_transcript, exons_n, length.
   - isoform_matrix.tsv: unique_jc, expression values per sample (or 1/0 if the expression files were not provided).
7. Create plots and summary tables (`sq_compare_summary.py`). Generate plots and tables summarizing isoform data, see /test/example_output.

---

## Output

/normalized_expression: a folder containing the normalized expression values if provided.
/summarized: a folder with the output tables and plots.

---

## Pipeline parameters

--input: tab-separated file with paths to output SQANTI3 files

--out: output directory

--collapseISM (optional)

Example:
`python sq_compare.py --input_files /path/to/sq_input_files.txt --out /path/to/output/folder'
