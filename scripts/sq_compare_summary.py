# report_isoforms.py

import pandas as pd
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib import rcParams
import seaborn as sns
from upsetplot import UpSet, from_memberships
import argparse

import warnings
import logging

matplotlib.use('Agg')  # Use non-interactive backend to avoid Qt errors
warnings.filterwarnings("ignore")
logging.getLogger('matplotlib').setLevel(logging.ERROR)
plt.ioff()

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
    "Fusion":"#FFB90F",  "Intergenic":"#E9967A", "Genic\nIntron":"#41B6C4"
}

grad_palette = ["#15918A", "#F58A53", "#FDC659"]  # green → orange → yellow

my_theme = {
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



parser = argparse.ArgumentParser(description="Summarize isoform metadata and expression matrices")
parser.add_argument(
    "--out", required=True,
    help="Output folder"
)
args = parser.parse_args()

isoform_info = pd.read_csv(f"{args.out}/summarized/isoform_info.tsv", sep="\t")
matrix = pd.read_csv(f"{args.out}/summarized/isoform_matrix.tsv", sep="\t")

samples = matrix.columns[1:].tolist() 
n_samples = len(samples)

# Summary statistics

summary_stats = {}

for sample in samples:
    present_jc = matrix.loc[matrix[sample] > 0, "unique_jc"]
    per_sample = isoform_info[isoform_info["unique_jc"].isin(present_jc)].copy()
    # total isoforms
    summary_stats.setdefault("isoforms_per_sample", pd.Series(dtype=int))
    summary_stats["isoforms_per_sample"].at[sample] = per_sample.shape[0]
    # isoforms per category
    cat_counts = per_sample["category"].value_counts().reindex(categories).fillna(0).astype(int)
    summary_stats.setdefault("isoforms_per_category", pd.DataFrame(0, index=samples, columns=categories))
    summary_stats["isoforms_per_category"].loc[sample] = cat_counts
    # mono vs multi
    mm_counts = per_sample["exons_n"].apply(lambda x: "Monoexon" if x==1 else "Multiexon").value_counts()
    summary_stats.setdefault("mono_vs_multi", pd.DataFrame(0, index=samples, columns=["Monoexon","Multiexon"]))
    summary_stats["mono_vs_multi"].loc[sample] = mm_counts

# Prepare summary tables for PDF
summary_rows = []

# Table 1: Isoforms per Sample
isoforms_per_sample = summary_stats["isoforms_per_sample"]
isoforms_per_sample_df = pd.DataFrame({
    "Sample": isoforms_per_sample.index,
    "Isoform count": isoforms_per_sample.values
})

# Table 2: Isoforms per Category
cat_table = summary_stats["isoforms_per_category"].copy()
cat_table.index.name = "Sample"
cat_table.reset_index(inplace=True)

# Table 3: Mono- VS Multiexon
mono_multi_table = summary_stats["mono_vs_multi"].copy()
mono_multi_table.index.name = "Sample"
mono_multi_table.reset_index(inplace=True)

# Write summary tables to PDF
report_path = f"{args.out}/summarized/sq_compare_report.pdf"
with PdfPages(report_path) as pdf:
    # Front page
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.7, "SQcompare Summary", fontsize=28, fontweight='bold', ha='center')
    ax.text(0.5, 0.6, f"Samples: {n_samples}", fontsize=18, ha='center')
    ax.text(0.5, 0.55, f"Total Isoforms: {isoform_info.shape[0]}", fontsize=16, ha='center')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    # Summary page
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.set_title("Table 1. Isoforms per Sample", fontsize=16, pad=40)
    plt.subplots_adjust(top=0.85) 
    table1 = ax.table(
        cellText=isoforms_per_sample_df.values,
        colLabels=isoforms_per_sample_df.columns,
        loc='center'
    )
    table1.auto_set_font_size(False)
    table1.set_fontsize(12)
    table1.scale(1.2, 1.2)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    #Table2
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.set_title("Table 2. Isoforms per Category", fontsize=16, pad=40)
    plt.subplots_adjust(top=0.85) 
    table2 = ax.table(
        cellText=cat_table.values,
        colLabels=cat_table.columns,
        loc='center',
        cellLoc='right'
    )
    table2.auto_set_font_size(False)
    table2.set_fontsize(12)
    table2.scale(1.2, 1.2)
    for (row, col), cell in table2.get_celld().items():
        if row == 0:  # header row
            cell.set_height(cell.get_height() * 1.8) 
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    #Table3
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.set_title("Table 3. Mono- VS Multiexon", fontsize=16, pad=40)
    plt.subplots_adjust(top=0.85) 
    table3 = ax.table(
        cellText=mono_multi_table.values,
        colLabels=mono_multi_table.columns,
        loc='center',
        cellLoc='right'
    )
    table3.auto_set_font_size(False)
    table3.set_fontsize(12)
    table3.scale(1.2, 1.2)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Generate plots

if n_samples < 7:
    cat_counts_df = summary_stats["isoforms_per_category"].copy()
    fig, ax = plt.subplots(figsize=(8, 6))
    cat_counts_df.plot(
        kind="bar",
        stacked=False,
        color=[cat_palette[c] for c in categories],
        ax=ax
    )
    ax.set_ylabel("UJC count")
    ax.set_xlabel("Sample")
    ax.set_title("Categories per Sample")
    ax.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # Save the figure as a JPEG file
    fig.savefig(f"{args.out}/summarized/ujc_per_category.jpeg", format="jpeg", dpi=300, bbox_inches='tight')

if n_samples < 7:
    # Stacked bar plot (proportions)
    cat_props_df = cat_counts_df.div(cat_counts_df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(8, 6))
    bottom = np.zeros(n_samples)
    for cat in categories:
        ax.bar(
            cat_props_df.index,
            cat_props_df[cat],
            bottom=bottom,
            color=cat_palette[cat],
            label=cat
        )
        bottom += cat_props_df[cat].values
    ax.set_ylabel("Proportion of UJC Categories")
    ax.set_xlabel("Sample")
    ax.set_title("UJC Category Composition Per Sample", fontsize=15)
    ax.set_ylim(0, 1)
    ax.legend(
        title="Category",
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0.
    )
    fig.tight_layout(rect=[0, 0, 0.8, 1])
    fig.savefig(f"{args.out}/summarized/ujc_category_composition.jpeg", format="jpeg", dpi=300, bbox_inches='tight')
    plt.close(fig)

# average length distribution
fig, ax = plt.subplots(figsize=(8,6))
for sample in samples:
    subset = isoform_info[isoform_info["unique_jc"].isin(
        matrix.loc[matrix[sample]>0, "unique_jc"]
    )]
    if not subset["average_length"].empty:
        sns.kdeplot(subset["average_length"], label=sample, ax=ax)
    ax.set_title("UJC Length Distributions")
    ax.set_xlabel("Length")
    ax.set_ylabel("Density")
    ax.legend(
        title="Sample",
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0.
    )
    fig.tight_layout(rect=[0, 0, 0.8, 1])
    fig.savefig(f"{args.out}/summarized/ujc_length_distribution.jpeg", format="jpeg", dpi=300, bbox_inches='tight')
    plt.close(fig)

    # heatmap (if expression values provided)
if (matrix[samples].values > 1).any():
    variances = matrix[samples].var(axis=1)
    top = matrix.loc[variances.nlargest(1000).index] 
    data = top.set_index("unique_jc")[samples]
    # Log-transform the expression values (add 1 to avoid log(0))
    data_log = np.log1p(data)
    # Use seaborn clustermap for hierarchical clustering
    cg = sns.clustermap(
        data_log,
        cmap="YlGnBu",
        figsize=(10, 6),
        cbar_kws={"label": "Log(TMM+1)"},
        yticklabels=False,  # Hide isoform IDs
        xticklabels=True    # Show sample names
    )
    # Titles and labels
    cg.ax_heatmap.set_title("Log-transformed Expression Clustermap (TMM) - Top 1000 Variable UJCs", pad=20)
    cg.ax_heatmap.set_xlabel("Samples")
    cg.ax_heatmap.set_ylabel("Isoforms")
    # Save figure
    cg.savefig(
        f"{args.out}/summarized/expression_clustermap.jpeg",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(cg.fig)  # close the figure to free memory

# UpSet plot with category-colored bars (if n_samples < 7)
if n_samples < 7:
    memberships = []
    jc_to_cat = isoform_info.set_index("unique_jc")["category"].map(cat_map).to_dict()
    jc_to_cat_raw = isoform_info.set_index("unique_jc")["category"].to_dict()
    for jc, row in matrix.set_index("unique_jc")[samples].iterrows():
        present = [s for s in samples if row[s] > 0]
        memberships.append(tuple(present))
    data = from_memberships(memberships)
    # Map each index (combination) to the list of unique_jc in that group
    index_to_jc = {}
    for idx, present in zip(data.index, data.index):
        # Find all unique_jc that match this membership
        mask = (matrix[samples] > 0)
        for i, s in enumerate(samples):
            mask[s] = mask[s] if present[i] else ~mask[s]
        mask_all = mask.all(axis=1)
        index_to_jc[present] = matrix.loc[mask_all, "unique_jc"].tolist()
    # For each bar, compute category proportions
    bar_cat_props = []
    for present in data.index:
        jcs = index_to_jc[present]
        cats = [jc_to_cat[jc] for jc in jcs if jc in jc_to_cat]
        if cats:
            cat_counts = pd.Series(cats).value_counts().reindex(categories, fill_value=0)
            props = cat_counts / cat_counts.sum()
        else:
            props = pd.Series([0]*len(categories), index=categories)
        bar_cat_props.append(props)
    # Plot upset
    fig = plt.figure(figsize=(8,6))
    upset = UpSet(data, subset_size='count', show_counts=True)
    upset.plot(fig=fig)
    # Color the bars by category proportions
    ax = upset.axs['intersections']
    for i, (bar, props) in enumerate(zip(ax.patches, bar_cat_props)):
        bottom = bar.get_y()
        left = bar.get_x()
        width = bar.get_width()
        height = bar.get_height()
        y0 = bottom
        for cat in categories:
            h = height * props[cat]
            if h > 0:
                rect = Rectangle((left, y0), width, h, color=cat_palette[cat], linewidth=0)
                ax.add_patch(rect)
                y0 += h
        # Hide the original bar
        bar.set_visible(False)
    ax.set_title("UpSet Plot (bars colored by category composition)")
    fig.tight_layout()
    fig.savefig(f"{args.out}/summarized/upset_category_colored.jpeg", format="jpeg", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Plot a standard UpSet plot (no category coloring)
if n_samples < 7:
    memberships = []
    for jc, row in matrix.set_index("unique_jc")[samples].iterrows():
        present = [s for s in samples if row[s] > 0]
        memberships.append(tuple(present))
    data = from_memberships(memberships)
    fig = plt.figure(figsize=(8,6))
    upset = UpSet(data, subset_size='count', show_counts=True)
    upset.plot(fig=fig)
    fig.tight_layout()
    fig.savefig(f"{args.out}/summarized/upset_standard.jpeg", format="jpeg", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Plot proportions of mono- and multi-exons per sample
mono_multi_props = mono_multi_table.iloc[:, 1:].div(mono_multi_table.iloc[:, 1:].sum(axis=1), axis=0)
mono_multi_props.index = mono_multi_table["Sample"] if "Sample" in mono_multi_table.columns else mono_multi_table.index
fig, ax = plt.subplots(figsize=(8, 6))
mono_multi_props.plot(
    kind="bar",
    stacked=True,
    ax=ax,
    color=["#FDC659", "#3B125A"]
)
ax.set_ylabel("Proportion of Isoforms")
ax.set_xlabel("Sample")
ax.set_title("Proportion of Monoexon vs Multiexon Isoforms per Sample")
ax.legend(title="Exon Type", bbox_to_anchor=(1.05, 1), loc='upper left')
fig.tight_layout(rect=[0, 0, 0.8, 1])
fig.savefig(f"{args.out}/summarized/mono_multi_proportion.jpeg", format="jpeg", dpi=300, bbox_inches='tight')
plt.close(fig)