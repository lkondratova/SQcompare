#!/usr/bin/env nextflow

/*
UJC comparison pipeline (Nextflow version of sq_compare.py)

conda activate sq_compare
nextflow run main.nf --input_files samples.tsv --outdir results

 */

nextflow.enable.dsl=2

params.input_files = file(params.input_files ?: "inputs.tsv")
params.out         = params.out ?: "results"
params.collapseISM = params.collapseISM ?: false

workflow {

    // 1: Parse inputs
    parse_out = parse(params.input_files, params.out)

    // 2: Collapse ISM (optional)
    if (params.collapseISM) {
        collapse_out = collapse_ism(parse_out.out_pickle, params.out)
        pickle_df    = collapse_out.out_pickle
    } else {
        pickle_df    = parse_out.out_pickle
    }

    // 3: Assign universal IDs
    universal_out = universal_id(pickle_df, params.out)

    // 4: Check expression file format and normalize if needed
    check_cols = count_columns(params.input_files)

    if (check_cols.num_cols == 4) {
        norm_out = tmm_norm(universal_out.out_pickle, params.out)
        matrix_out = generalize_isoforms(norm_out.out_pickle, params.out)
    } else {
        matrix_out = generalize_isoforms(universal_out.out_pickle, params.out)
    }

    // 5: Visualization / summary
    compare_out = sq_compare_summary(params.out)

    // 6: Export report
    export_out = export_report(params.out)
}


// --------------------------- Processes ---------------------------

// 1. Parse SQANTI3 inputs
process parse {
    input:
    path input_files
    val outdir

    output:
    path "${outdir}/sqanti3_samples.pkl" into parsed_pickle

    script:
    """
    python scripts/parse_sq_inputs.py --input_files ${input_files} --out ${outdir}
    """
}

// 2. Collapse ISM isoforms
process collapse_ism {
    input:
    path pickle_file
    val outdir

    output:
    path "${outdir}/sqanti3_samples_ISMcollapsed.pkl" into collapsed_pickle

    script:
    """
    python scripts/collapse_ism.py --pickle ${pickle_file} --out ${outdir}
    """
}

// 3. Assign universal IDs
process universal_id {
    input:
    path pickle_file
    val outdir

    output:
    path "${outdir}/sqanti3_standardized.pkl" into standardized_pickle

    script:
    """
    python scripts/universal_id.py --pickle ${pickle_file} --out ${outdir}
    """
}

// 4. Count columns in input file to decide if expression values are present
process count_columns {
    input:
    path input_files

    output:
    tuple val(input_files) val(num_cols)

    script:
    """
    NUM_COLS=\$(head -n 1 ${input_files} | awk -F'\\t' '{print NF}')
    echo "${input_files}\t\${NUM_COLS}" > cols.txt
    """
}

// 5a. Normalize expression values with edgeR TMM
process tmm_norm {
    input:
    path pickle_file
    val outdir

    output:
    path "${outdir}/sqanti3_normalized.pkl" into normalized_pickle

    script:
    """
    python scripts/tmm_norm.py --pickle ${pickle_file} --out ${outdir}
    """
}

// 5b. Create isoform matrix and info files
process generalize_isoforms {
    input:
    path pickle_file
    val outdir

    output:
    path "${outdir}/isoform_matrix.tsv"
    path "${outdir}/isoform_info.tsv"
    path "${outdir}/summary_stats.tsv"

    script:
    """
    python scripts/generalize_isoforms.py --pickle ${pickle_file} --out ${outdir}
    """
}

// 6. Generate summary plots
process sq_compare_summary {
    input:
    val outdir

    output:
    path "${outdir}/plots/"

    script:
    """
    python scripts/sq_compare_summary.py --out ${outdir}
    """
}

// 7. Export final report (PDF/TSV/whatever export_script produces)
process export_report {
    input:
    val outdir

    output:
    path "${outdir}/final_report.pdf"

    script:
    """
    python scripts/export_script.py --out ${outdir}
    """
}
