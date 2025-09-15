# SQcompare

This repository contains a Python-based workflow to analyze and compare multiple SQANTI3 outputs. The pipeline handles multiple samples, collapses incomplete splice match (ISM) isoforms (optional), assigns universal isoform IDs across samples, normalizes expression, and generates summary plots and reports.

---

## Features

- Parse SQANTI3 output files: classification, junctions, GTF, expression (optional).
- Collapse ISM isoforms (optional).
- Assign **universal isoform IDs** across samples based on junction chains.
- Normalize expression values using TMM (edgeR-like) normalization if expression files are provided.
- Generate combined isoform matrices and summary files.
- Produce publication-ready plots and tables:
  - Isoform counts per category
  - Length distributions
  - Heatmaps of expression
  - UpSet plots for isoform sharing
  - Monoexon vs multiexon counts
- Export a final **PDF report** containing plots, tables, and summary statistics.

---

## Installation

This pipeline is designed to run in a **Conda environment**. Use the provided `sq_compare_environment.yml`:

```bash
conda env create -f sq_compare_environment.yml
conda activate sq_compare
```

Nextflow to automate the pipeline
`nextflow run main.nf -profile standard`

---
## Input Files

SQcompare accepts a tab-separated input file (no header) with the following columns:

1. Path to a SQANTI3 *classification.txt

2. Path to a SQANTI3 *junctions.txt

3. Path to *corrected.gtf

4. Optionally, an expression file (absolute values) can be provided, where the first column is an isoform name and the second column is an absolute expression value (no header).

An example bash script for generating an input file can be found in `helper_scripts/create_input_file.sh`

---

## Workflow Overview (Nextflow)

1. Parse inputs: `parse_sq_inputs.py` (Parse SQANTI3 files and organize them into dataframes)
2. Collapse ISM isoforms (optional): `collapse_ism.py` (Collapse incomplete splice match isoforms if requested)
3. Assign universal IDs: `universal_id.py` (Assign universal isoform IDs across all samples based on junction chains)
4. Normalize expression (optional): `tmm_norm.py` (Normalize expression values using TMM edgeR-like normalization)
5. Generate matrices and combined isoform info: `generalize_isoforms.py` (Create combined isoform matrices and information files for all samples)
6. Create plots and summary tables: `sq_compare_summary.py` (Generate plots and tables summarizing isoform data)
7. Export report: `export_script.py` (Export plots, tables, and summary statistics, optionally as a single PDF report)

---

## Output

`isoform_info.tsv`: Summary of all isoforms with universal IDs, category, length, exons, etc.

`isoform_matrix.tsv`: Combined matrix with isoform presence/absence or normalized expression per sample.

`plots/`: Figures for category counts, length distributions, heatmaps, UpSet plots, and exon distribution.

`tables/`: TSV files with counts and summary tables.

`summary_report.txt`: Text summary of isoform statistics.

`full_report.pdf`: Optional combined report with all plots, tables, and statistics.

---

## Configuration

Nextflow config: nextflow.config

Profiles for local, cluster, or Docker execution.

Uses the provided sq_compare_environment.yml for reproducibility.

## Pipeline parameters

--input: tab-separated file with paths to output SQANTI3 files

--out: output directory

--collapseISM (optional)

## Dependencies

Python ≥ 3.10

pandas, numpy, scipy

matplotlib, seaborn

upsetplot

reportlab

Conda (recommended for reproducibility)
