# export_script.py
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from report_script import (
    plot1_bar, plot1_table,
    plot2_dist,
    plot3_heatmap,
    plot4_upset,
    plot5_bar, plot5_table,
    summary_stats
)

def export_report(output_dir="report_output", make_pdf=True):
    """
    Export all plots, tables, and summary statistics.
    If make_pdf=True, also export a combined PDF report.
    """

    os.makedirs(output_dir, exist_ok=True)

    pdf_path = os.path.join(output_dir, "full_report.pdf")
    pdf = PdfPages(pdf_path) if make_pdf else None

    # === Plots ===
    def save_plot(plot, name):
        if plot is not None:
            fig = plot.figure if hasattr(plot, "figure") else plot
            fig.savefig(os.path.join(output_dir, f"{name}.png"), dpi=300, bbox_inches="tight")
            if pdf:
                pdf.savefig(fig)

    save_plot(plot1_bar, "plot1_bar")
    save_plot(plot2_dist, "plot2_length_distribution")
    save_plot(plot3_heatmap, "plot3_heatmap")
    save_plot(plot4_upset, "plot4_upset")
    save_plot(plot5_bar, "plot5_bar")

    # === Tables ===
    def save_table(table, name):
        if table is not None:
            table.to_csv(os.path.join(output_dir, f"{name}.tsv"), sep="\t", index=True)
            if pdf:
                fig, ax = plt.subplots(figsize=(8, len(table) * 0.3 + 1))
                ax.axis("off")
                ax.table(
                    cellText=table.reset_index().values,
                    colLabels=table.reset_index().columns,
                    loc="center"
                )
                ax.set_title(name, fontsize=12, fontweight="bold")
                pdf.savefig(fig)
                plt.close(fig)

    save_table(plot1_table, "table1_isoforms_per_category")
    save_table(plot5_table, "table5_isoforms_per_sample")

    # === Summary statistics ===
    if summary_stats is not None:
        # Save TXT file
        with open(os.path.join(output_dir, "summary_report.txt"), "w") as f:
            f.write("=== Isoform Report Summary ===\n\n")
            for key, val in summary_stats.items():
                f.write(f"{key}:\n{val}\n\n")

        if pdf:
            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")
            text = "=== Isoform Report Summary ===\n\n"
            for key, val in summary_stats.items():
                text += f"{key}:\n{val}\n\n"
            ax.text(0.01, 0.99, text, va="top", ha="left", fontsize=10, family="monospace")
            pdf.savefig(fig)
            plt.close(fig)

    # === Close PDF if created ===
    if pdf:
        pdf.close()
        print(f"PDF report saved at: {pdf_path}")

    print(f"Report exported to: {output_dir}")


if __name__ == "__main__":
    export_report()
