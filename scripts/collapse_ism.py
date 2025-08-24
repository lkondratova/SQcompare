import os
import argparse
import pickle
import pandas as pd
from collections import defaultdict

def collapse_ISM(class_df, junc_df, expr_df=None):
    collapsed_dict = {}
    kept = []
    dropped = []

    # Group by associated_transcript
    for at, group in class_df.groupby("associated_transcript"):
        if at == "novel":
            kept.extend(group["isoform"].tolist())
            continue

        if len(group) == 1:
            kept.extend(group["isoform"].tolist())
            continue

        # case handling
        sc = group["structural_category"].tolist()
        isoforms = group["isoform"].tolist()

        # Case 1: FSM present
        if "full_splice_match" in sc:
            fsm = group[group["structural_category"] == "full_splice_match"]

            if len(fsm) == 1:
                survivor = fsm.iloc[0]["isoform"]
            else:
                # multiple FSMs -> use subcategory priority
                priority = ["reference_match", "alternative_5end", "alternative_3end", "alternative_3end5end"]
                fsm_sorted = fsm.sort_values(
                    by="subcategory",
                    key=lambda col: col.map(lambda x: priority.index(x) if x in priority else len(priority))
                )
                survivor = fsm_sorted.iloc[0]["isoform"]

            removed = [i for i in isoforms if i != survivor]
            collapsed_dict[survivor] = removed
            kept.append(survivor)
            dropped.extend(removed)

        # Case 2: Only ISMs
        else:
            priority = ["5prime_fragment", "3prime_fragment", "internal_fragment"]
            ism_sorted = group.sort_values(
                by="subcategory",
                key=lambda col: col.map(lambda x: priority.index(x) if x in priority else len(priority))
            )
            survivor = ism_sorted.iloc[0]["isoform"]
            removed = [i for i in isoforms if i != survivor]
            collapsed_dict[survivor] = removed
            kept.append(survivor)
            dropped.extend(removed)

    # Update classification
    class_df = class_df[class_df["isoform"].isin(kept)]

    # Update junctions
    junc_df = junc_df[~junc_df["isoform"].isin(dropped)]

    # Update expression
    if expr_df is not None:
        expr_df = expr_df.copy()
        iso_to_expr = dict(zip(expr_df.iloc[:, 0], expr_df.iloc[:, 1]))

        for survivor, removed in collapsed_dict.items():
            expr_df.loc[expr_df.iloc[:, 0] == survivor, expr_df.columns[1]] += sum(iso_to_expr.get(r, 0) for r in removed)

        expr_df = expr_df[~expr_df.iloc[:, 0].isin(dropped)]

    return class_df, junc_df, expr_df, collapsed_dict


def main():
    parser = argparse.ArgumentParser(description="Collapse ISM isoforms per sample")
    parser.add_argument("--class", required=True, help="Classification file")
    parser.add_argument("--junctions", required=True, help="Junctions file")
    parser.add_argument("--expression", help="Expression file (optional)")
    parser.add_argument("--pickle", required=True, help="Pickle file containing parsed dataframes")
    parser.add_argument("--out", required=True, help="Output folder")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load files
    class_df = pd.read_csv(args.__dict__["class"], sep="\t")
    junc_df = pd.read_csv(args.junctions, sep="\t")
    expr_df = pd.read_csv(args.expression, sep="\t") if args.expression else None

    # Collapse
    class_df, junc_df, expr_df, collapsed_dict = collapse_ISM(class_df, junc_df, expr_df)

    # Save new files
    def save_file(df, in_path, suffix="ISMcollapsed"):
        name, ext = os.path.splitext(os.path.basename(in_path))
        out_path = os.path.join(args.out, f"{name}.{suffix}{ext}")
        df.to_csv(out_path, sep="\t", index=False)
        return out_path

    class_out = save_file(class_df, args.__dict__["class"])
    junc_out = save_file(junc_df, args.junctions)
    expr_out = save_file(expr_df, args.expression) if expr_df is not None else None

    # Save collapsed summary
    summary_path = os.path.join(args.out, "ISMcollapsed_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("survivor_isoform\tcollapsed_isoforms\n")
        for survivor, removed in collapsed_dict.items():
            if removed:  # only log when something collapsed
                f.write(f"{survivor}\t{','.join(removed)}\n")

    # Update pickle so downstream scripts use collapsed data
    with open(args.pickle, "rb") as f:
        data = pickle.load(f)

    data["classification"] = class_df
    data["junctions"] = junc_df
    if expr_df is not None:
        data["expression"] = expr_df

    with open(args.pickle, "wb") as f:
        pickle.dump(data, f)

    print(f"Collapsed ISMs. Updated files saved to {args.out}. Pickle updated.")
    print(f"Collapsed summary written to {summary_path}")


if __name__ == "__main__":
    main()