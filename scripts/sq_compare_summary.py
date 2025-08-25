# report_isoforms.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from upsetplot import UpSet, from_memberships

# -------------------
# Settings
# -------------------

categories_raw = [
    "full-splice_match", "incomplete-splice_match",
    "novel_in_catalog", "novel_not_in_catalog",
    "genic", "antisense", "fusion", "intergenic", "genic_intron"
]

categories = [
    "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",
    "Antisense", "Fusion", "Intergenic", "Genic\nIntron"
]

cat_map = dict(zip(categories_raw, categories))

cat_palette = {
    "FSM":"#6BAED6", "ISM":"#FC8D59", "NIC":"#78C679", 
    "NNC":"#EE6A50", "Genic\nGenomic":"#969696", "Antisense":"#66C2A4", 
    "Fusion":"goldenrod1",  "Intergenic":"darksalmon", "Genic\nIntron":"#41B6C4"
}

grad_palette = ["#15918A", "#F58A53", "#FDC659"]  # green → orange → yellow

my_theme = {
    "font.family": "Helvetica",
    "axes.titlesize": 15,
    "axes.titleweight": "normal",
    "axes.titlepad": 10,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.title_fontsize": 12,
    "legend.fontsize": 11,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.4,
    "legend.handlelength": 1,   
    "legend.handleheight": 0.5, 
    "figure.subplot.left": 0.1,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.1,
    "figure.subplot.top": 0.9,
}
rcParams.update(my_theme)

# -------------------
# Load data
# -------------------

isoform_info = pd.read_csv("isoform_info.tsv", sep="\t")
matrix = pd.read_csv("isoform_matrix.tsv", sep="\t")

samples = matrix.columns[1:].tolist()  # assuming first column is unique_jc
n_samples = len(samples)

# -------------------
# Summary statistics
# -------------------

summary_stats = {}

# total isoforms
summary_stats["total_isoforms"] = isoform_info.shape[0]

# isoforms per category
summary_stats["isoforms_per_category"] = (
    isoform_info["category"].map(cat_map).value_counts()[categories]
)

# isoforms per sample
summary_stats["isoforms_per_sample"] = matrix[samples].astype(bool).sum()

# average isoform length
summary_stats["average_length"] = isoform_info["average_length"].mean()

# monoexon vs multiexon
summary_stats["mono_vs_multi"] = isoform_info["exons_n"].apply(lambda x: "Monoexon" if x==1 else "Multiexon").value_counts()

# expression summary if available
if (matrix[samples].values > 1).any():  # heuristic: expression values provided
    summary_stats["expression_range"] = {
        "min": matrix[samples].min().min(),
        "max": matrix[samples].max().max(),
        "mean": matrix[samples].mean().mean()
    }

# -------------------
# Plots & Tables
# -------------------

plots = {}
tables = {}

# ---- plot1: barplot isoforms per category/sample ----
cat_table = pd.DataFrame(0, index=samples, columns=categories)
for sample in samples:
    subset = isoform_info[isoform_info["unicue_jc"].isin(
        matrix.loc[matrix[sample]>0, "unicue_jc"]
    )]
    counts = subset["category"].map(cat_map).value_counts()
    cat_table.loc[sample, counts.index] = counts.values

tables["cat_table"] = cat_table

if n_samples <= 8:
    fig, ax = plt.subplots(figsize=(8,6))
    cat_table.plot(kind="bar", stacked=False, color=[cat_palette[c] for c in categories], ax=ax)
    ax.set_ylabel("Isoform count")
    ax.set_title("Isoforms per category per sample")
    plots["plot1"] = fig

# ---- plot2: average length distribution ----
fig, ax = plt.subplots(figsize=(8,6))
for sample in samples:
    subset = isoform_info[isoform_info["unicue_jc"].isin(
        matrix.loc[matrix[sample]>0, "unicue_jc"]
    )]
    sns.kdeplot(subset["average_length"], label=sample, ax=ax)
ax.set_title("Isoform length distributions")
ax.set_xlabel("Length")
ax.set_ylabel("Density")
ax.legend()
plots["plot2"] = fig

# ---- plot3: heatmap (if expression values provided) ----
if (matrix[samples].values > 1).any():
    fig, ax = plt.subplots(figsize=(10,6))
    sns.heatmap(matrix.set_index("unicue_jc")[samples], cmap="YlGnBu", ax=ax)
    ax.set_title("Expression heatmap (TMM)")
    plots["plot3"] = fig

# ---- plot4: UpSet (skip if too many samples) ----
if n_samples <= 6:
    memberships = []
    jc_to_cat = isoform_info.set_index("unicue_jc")["category"].map(cat_map).to_dict()
    for jc, row in matrix.set_index("unicue_jc")[samples].iterrows():
        present = [s for s in samples if row[s] > 0]
        memberships.append(present)
    data = from_memberships(memberships)
    fig = plt.figure(figsize=(8,6))
    upset = UpSet(data, subset_size='count', show_counts=True)
    upset.plot(fig=fig)
    plots["plot4"] = fig

# ---- plot5: mono vs multi per sample ----
mono_multi_table = pd.DataFrame(0, index=samples, columns=["Monoexon","Multiexon"])
for sample in samples:
    subset = isoform_info[isoform_info["unicue_jc"].isin(
        matrix.loc[matrix[sample]>0, "unicue_jc"]
    )]
    counts = subset["exons_n"].apply(lambda x: "Monoexon" if x==1 else "Multiexon").value_counts()
    mono_multi_table.loc[sample, counts.index] = counts.values

tables["mono_multi_table"] = mono_multi_table

if n_samples <= 8:
    fig, ax = plt.subplots(figsize=(8,6))
    mono_multi_table.plot(kind="bar", stacked=True, ax=ax, color=["#41B6C4","#EE6A50"])
    ax.set_ylabel("Isoform count")
    ax.set_title("Monoexon vs Multiexon isoforms per sample")
    plots["plot5"] = fig

# -------------------
# Outputs
# -------------------

report = {
    "summary_stats": summary_stats,
    "plots": plots,
    "tables": tables
}