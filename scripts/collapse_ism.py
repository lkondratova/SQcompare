import os
import argparse
import pickle
import pandas as pd
#from collections import defaultdict


def collapse_ISM(sample, class_df, junc_df, expr_df=None):
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
        expr_df['survivor'] = expr_df.iloc[:, 0].map(lambda x: next((s for s, r in collapsed_dict.items() if x in r), x))
        expr_df = expr_df.groupby('survivor').agg({expr_df.columns[0]: 'first', expr_df.columns[1]: 'sum'}).reset_index(drop=True)
        expr_df.rename(columns={'survivor': expr_df.columns[0]}, inplace=True)
        expr_df = expr_df[expr_df[expr_df.columns[0]].isin(kept)]
    sample_dict={sample: collapsed_dict}
    return class_df, junc_df, expr_df, sample_dict

# Save new files
def save_file(df, out, sample, suffix):
    out_path = os.path.join(out, f"{sample}_{suffix}.ISMcollapsed.txt")
    df.to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Collapse ISM isoforms per sample")
    parser.add_argument("--pickle", required=True, help="Pickle file containing parsed dataframes")
    parser.add_argument("--out", required=True, help="Output folder")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    with open(args.pickle, "rb") as f:
        data = pickle.load(f)
        for sample in data['samples']:
            class_df=data['data'][sample]['classification']
            junc_df=data['data'][sample]['junctions']
            gtf_df=data['data'][sample]['gtf']
            expr_df=data['data'][sample]['expression'] if 'expression' in data['data'][sample] else None
            
            # Collapse
            class_df, junc_df, expr_df, collapsed_dict = collapse_ISM(sample, class_df, junc_df, expr_df)
            #save_file(class_df, args.out, sample, "classification")
            #save_file(junc_df, args.out, sample, "junctions")
            #save_file(expr_df, args.out, sample, "expression") if expr_df is not None else None

            # Filter GTF
            attr=gtf_df[8]
            gtf_df['transcript_id'] = attr.str.extract('transcript_id "([^"]+)"')
            gtf_df = gtf_df[gtf_df['transcript_id'].isin(class_df['isoform'])]

            # Update sample data in pickle with collapsed files
            data['data'][sample]['classification'] = class_df
            data['data'][sample]['junctions'] = junc_df
            data['data'][sample]['gtf'] = gtf_df.drop(columns=['transcript_id'])
            if expr_df is not None:
                data['data'][sample]['expression'] = expr_df

            print(f"Processed sample {sample}: Kept {len(class_df)} isoforms.")  

    # Save collapsed summary
    summary_path = os.path.join(args.out, "ISMcollapsed_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\tsurvivor_isoform\tcollapsed_isoforms\n")
        for sample, transcripts in collapsed_dict.items():
            for survivor, removed in transcripts.items():
                if removed:  # only log when something collapsed
                    f.write(f"{sample}\t{survivor}\t{','.join(removed)}\n")
   
    # Replace old pickle 
    with open(f'{args.out}/sqanti3_samples_ISMcollapsed.pkl', "wb") as f:
        pickle.dump(data, f)

    #print(f"Collapsed ISMs. Updated files saved to {args.out}. Pickle updated.")
    print(f"Collapsed summary written to {summary_path}")


if __name__ == "__main__":
    main()