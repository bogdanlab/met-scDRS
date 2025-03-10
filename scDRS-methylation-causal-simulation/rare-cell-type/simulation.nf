#!/usr/bin/env nextflow
/*
 * pipeline for simulate causal perturbation to investigate power for different cell types
 */

// generate cell type meta data:
process get_cell_type_metadata {
    input:
    //input files:
    path executable
    val meta_path
    val output_path

    output:
    // output file:
    val output_path

    script:
    """
    Rscript ${executable} \
        --meta_csv ${meta_path} \
        --meta_column "MajorType" \
        --output ${output_path}
    """
}

// define perturbation processes that generate simulated perturbed expression with set overlap:
process perturb_fraction_VLMC {
    input:
    // input files:
    path executable
    path subset_file
    path cell_type_meta
    path gs_file
    val cell_type
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --cell_type_meta ${cell_type_meta} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_type ${cell_type} \
            --proportion_to_perturb 0.25 \
            --down_sample_to 0.2 \
            --down_sample_rate 0.2 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 10 \
            --output_dir ${output_dir}
    """
}
process perturb_fraction_L6b {
    input:
    // input files:
    val pass
    path executable
    path subset_file
    path cell_type_meta
    path gs_file
    val cell_type
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --cell_type_meta ${cell_type_meta} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_type ${cell_type} \
            --proportion_to_perturb 0.25 \
            --down_sample_to 0.2 \
            --down_sample_rate 0.2 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 10 \
            --output_dir ${output_dir}
    """
}
process perturb_fraction_ASC {
    input:
    // input files:
    val pass
    path executable
    path subset_file
    path cell_type_meta
    path gs_file
    val cell_type
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --cell_type_meta ${cell_type_meta} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_type ${cell_type} \
            --proportion_to_perturb 0.25 \
            --down_sample_to 0.2 \
            --down_sample_rate 0.2 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 10 \
            --output_dir ${output_dir}
    """
}
process perturb_fraction_ITL4 {
    input:
    // input files:
    val pass
    path executable
    path subset_file
    path cell_type_meta
    path gs_file
    val cell_type
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --cell_type_meta ${cell_type_meta} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_type ${cell_type} \
            --proportion_to_perturb 0.25 \
            --down_sample_to 0.2 \
            --down_sample_rate 0.2 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 10 \
            --output_dir ${output_dir}
    """
}
process perturb_fraction_ITL23 {
    input:
    // input files:
    val pass
    path executable
    path subset_file
    path cell_type_meta
    path gs_file
    val cell_type
    // directories for outputting data:
    val output_dir

    output:
    val output_dir

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --cell_type_meta ${cell_type_meta} \
            --gs_file ${gs_file} \
            --trait "UKB_460K.body_HEIGHTz" \
            --perturbation_effect_start 1.005 \
            --perturbation_effect_end 1.005 \
            --effect_step 0.01 \
            --gene_number 1000 \
            --cell_type ${cell_type} \
            --proportion_to_perturb 0.25 \
            --down_sample_to 0.2 \
            --down_sample_rate 0.2 \
            --overlap_start 0.25 \
            --overlap_end 0.25 \
            --overlap_step 0.1 \
            --replication 10 \
            --output_dir ${output_dir}
    """
}
// define parameters:
// for extracting meta information:
params.meta_file = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
params.get_meta_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-causal-simulation/rare-cell-type/get_cell_type_meta_data.R"
params.output_meta_file = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/cell_type_meta_only.csv"

// for perturbation
// cell_type_perturb = ["VLMC", "L6b", "ASC", "IT-L4", "IT-L23"]

// for execute:
params.subsetted_file = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv"
params.gs_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
params.overlap_output_dir = "/u/scratch/l/lixinzhe/tmp-file/causal-simulation/rare-cell-type/"
params.executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-causal-simulation/rare-cell-type/perturbation.R"
params.cell_type_VLMC = "VLMC"
params.cell_type_L6b = "L6b"
params.cell_type_ASC = "ASC"
params.cell_type_ITL4 = "IT-L4"
params.cell_type_ITL23 = "IT-L23"

// define execution workflow:
workflow {
    meta_file = get_cell_type_metadata(params.get_meta_executable, params.meta_file, params.output_meta_file)
    perturbed_VLMC = perturb_fraction_VLMC(params.executable, params.subsetted_file, meta_file, params.gs_file, params.cell_type_VLMC, params.overlap_output_dir)
    perturbed_L6b = perturb_fraction_L6b(perturbed_VLMC, params.executable, params.subsetted_file, meta_file, params.gs_file, params.cell_type_L6b, params.overlap_output_dir)
    perturbed_ASC = perturb_fraction_ASC(perturbed_L6b, params.executable, params.subsetted_file, meta_file, params.gs_file, params.cell_type_ASC, params.overlap_output_dir)
    perturbed_ITL4 = perturb_fraction_ITL4(perturbed_ASC, params.executable, params.subsetted_file, meta_file, params.gs_file, params.cell_type_ITL4, params.overlap_output_dir)
    perturbed_ITL23 = perturb_fraction_ITL23(perturbed_ITL4, params.executable, params.subsetted_file, meta_file, params.gs_file, params.cell_type_ITL23, params.overlap_output_dir)
}
