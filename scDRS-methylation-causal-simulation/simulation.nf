#!/usr/bin/env nextflow
/*
 * pipeline for methylation causal simulation:
 */

// define getting the perturbed expression matrix

// define perturbation processes that generate simulated perturbed expression with set overlap:
process causal_perturbation_fix_overlap {
    input:
    // input files:
    path executable
    path subset_file
    path gs_file
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.01 \
            --perturbation_effect_end 1.05 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_number 500 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 100 \
            --output_dir ${output_dir}
    """
}

process causal_perturbation_fix_overlap_small {
    input:
    // input files:
    path executable
    path subset_file
    path gs_file
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.001 \
            --perturbation_effect_end 1.01 \
            --effect_step 0.001 \
            --gene_number 1000 \
            --cell_number 500 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 100 \
            --output_dir ${output_dir}
    """
}


// define perturbation processes that generate simulated
process causal_perturbation_fix_effect {
    input:
    // input files:
    path executable
    path subset_file
    path gs_file
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.001 \
            --gene_number 1000 \
            --cell_number 500 \
            --overlap_start 0.1 \
            --overlap_end 0.5 \
            --overlap_step 0.1 \
            --replication 100 \
            --output_dir ${output_dir}
    """
}

// for perturbation::
params.subsetted_file = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv"
params.gs_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
params.overlap_output_dir = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-overlap/"
params.effect_output_dir = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-effect/"
params.executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-causal-simulation/perturbation.R"

// define execution workflow:
workflow {
    //fix_overlap_perturbed_small_dir = causal_perturbation_fix_overlap_small(params.executable, params.subsetted_file, params.gs_file, params.overlap_output_dir)
    //fix_overlap_perturbed_dir = causal_perturbation_fix_overlap(fix_overlap_perturbed_small_dir, params.executable, params.subsetted_file, params.gs_file, params.overlap_output_dir)
    fix_effect_perturbed_dir = causal_perturbation_fix_effect(params.executable, params.subsetted_file, params.gs_file, params.effect_output_dir)
}
