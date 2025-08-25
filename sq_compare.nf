#!/usr/bin/env nextflow

/*
 * Isoform analysis workflow
 */

params.input           = "./input/*"        // path to raw input files
params.samples         = "./samples.txt"    // optional sample annotation
params.collapseISM     = true               // collapse ISM or not
params.expression      = "./expression/*"   // path to expression files (optional)
params.outdir          = "./results"        // output folder

process PARSE {
    tag "$sample_id"

    input:
    path input_files from channel.fromPath(params.input)

    output:
    path "parsed/" into parsed_ch

    script:
    """
    mkdir -p parsed
    python parse_sq_inputs.py --input $input_files --output parsed
    """
}

process COLLAPSE_ISM {
    when:
    params.collapseISM

    input:
    path parsed_files from parsed_ch

    output:
    path "collapsed/" into collapsed_ch

    script:
    """
    mkdir -p collapsed
    python collapse_ism.py --input $parsed_files --output collapsed
    """
}

process UNIVERSAL_ID {
    input:
    path infiles from (params.collapseISM ? collapsed_ch : parsed_ch)

    output:
    path "universal/" into uid_ch

    script:
    """
    mkdir -p universal
    python universal_id.py --input $infiles --output universal
    """
}

process TMM_NORM {
    when:
    params.expression

    input:
    path expr_files from channel.fromPath(params.expression)
    path isoform_files from uid_ch

    output:
    path "normalized/" into norm_ch

    script:
    """
    mkdir -p normalized
    python tmm_norm.py --input $expr_files --isoforms $isoform_files --output normalized
    """
}

process GENERALIZE {
    input:
    path infiles from (params.expression ? norm_ch : uid_ch)

    output:
    path "generalized/" into gen_ch

    script:
    """
    mkdir -p generalized
    python generalize_isoforms.py --input $infiles --output generalized
    """
}

process PLOTS {
    input:
    path gen_files from gen_ch

    output:
    path "plots/" into plots_ch

    script:
    """
    mkdir -p plots
    python sq_compare_summary.py --input $gen_files --output plots
    """
}

process EXPORT {
    input:
    path plots from plots_ch

    output:
    path "${params.outdir}/final_report"

    script:
    """
    mkdir -p ${params.outdir}/final_report
    python export_script.py --input $plots --output ${params.outdir}/final_report
    """
}

workflow {
    PARSE()
    if (params.collapseISM) {
        COLLAPSE_ISM(PARSE.out)
        UNIVERSAL_ID(COLLAPSE_ISM.out)
    } else {
        UNIVERSAL_ID(PARSE.out)
    }

    if (params.expression) {
        TMM_NORM(UNIVERSAL_ID.out)
        GENERALIZE(TMM_NORM.out)
    } else {
        GENERALIZE(UNIVERSAL_ID.out)
    }

    PLOTS(GENERALIZE.out)
    EXPORT(PLOTS.out)
}
