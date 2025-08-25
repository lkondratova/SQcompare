# SQcompare

This repository contains a Python-based workflow to analyze and compare multiple SQANTI3 outputs. The pipeline handles multiple samples, collapses incomplete splice match (ISM) isoforms (optional), assigns universal isoform IDs across samples, normalizes expression, and generates summary plots and reports.

---

## Features

- Parse SQANTI3 output files (classification, junctions, GTF, optional expression).
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

This pipeline is designed to run in a **Conda environment**. Use the provided `environment.yml`:

```bash
conda env create -f environment.yml
conda activate isoform-tool

Nextflow to automate the pipeline
nextflow run main.nf -profile standard
Input Files

SQANTI3 output files:
Each sample should have at least:

classification.txt

junctions.txt

isoforms.gtf

Optionally, an expression file (TPM or counts) can be provided.

Samples TSV: (optional)
Can include metadata like sample names, days, or conditions.

Pipeline Scripts
Script	Purpose
parse_sq_inputs.py	Parse SQANTI3 files and organize them into dataframes.
collapse_ism.py	Collapse incomplete splice match isoforms if requested.
universal_id.py	Assign universal isoform IDs across all samples based on junction chains.
tmm_norm.py	Normalize expression values using TMM (edgeR-like normalization).
generalize_isoforms.py	Create combined isoform matrices and information files for all samples.
sq_compare_summary.py	Generate plots and tables summarizing isoform data.
export_script.py	Export plots, tables, and summary statistics, optionally as a single PDF report.
Workflow Overview (Nextflow)

Parse inputs: parse_sq_inputs.py

Collapse ISM isoforms (optional): collapse_ism.py

Assign universal IDs: universal_id.py

Normalize expression (optional): tmm_norm.py

Generate matrices and combined isoform info: generalize_isoforms.py

Create plots and summary tables: sq_compare_summary.py

Export report: export_script.py

Output

isoform_info.tsv: Summary of all isoforms with universal IDs, category, length, exons, etc.

isoform_matrix.tsv: Combined matrix with isoform presence/absence or normalized expression per sample.

plots/: Figures for category counts, length distributions, heatmaps, UpSet plots, and exon distribution.

tables/: TSV files with counts and summary tables.

summary_report.txt: Text summary of isoform statistics.

full_report.pdf: Optional combined report with all plots, tables, and statistics.

Configuration

Nextflow config: nextflow.config

Profiles for local, cluster, or docker execution.

Uses the provided environment.yml for reproducibility.

Pipeline parameters:

--input: path to input SQANTI3 files

--samples: optional metadata file

--collapseISM: true/false

--expression: path to expression files (optional)

--outdir: output directory

Dependencies

Python ≥ 3.10

pandas, numpy, scipy

matplotlib, seaborn

upsetplot

reportlab

Nextflow (optional, for pipeline automation)

Conda (recommended for reproducibility)

Future Extensions

Differential isoform expression (DIE) or differential transcript usage (DTU) analysis.

Functional impact analysis of isoform switching (protein domain changes, NMD predictions).

Integration with external tools like IsoformSwitchAnalyzeR for deeper splicing analysis.

License

MIT License

Contact

For questions, please contact [Your Name] or open an issue on GitHub.


