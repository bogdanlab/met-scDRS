#!/usr/bin/env nextflow
/*
 * pipeline for methylation data processing:
 */

// define getting the high fraction gene set:
process high_fraction_genes {
    input:
    path executable
    path subset_file
    path gs_file
    each gene_number
    val output_dir

    output:
    val 'pass'

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --gs_file ${gs_file} \
            --gene_number ${gene_number} \
            --quantile 0.95 \
            --replication 100 \
            --output_dir ${output_dir} > "gene_num_${gene_number}.logfile"

    """
}

// define process for getting high variance gene set:
process high_variance_genes {
    input:
    val pass
    path executable
    path subset_file
    path gs_file
    each gene_number
    val output_dir

    output:
    val 'pass'

    script:
    """
    Rscript ${executable} \
            --data_matrix ${subset_file} \
            --gs_file ${gs_file} \
            --gene_number ${gene_number} \
            --quantile 0.95 \
            --replication 100 \
            --output_dir ${output_dir} > "gene_num_${gene_number}.logfile"

    """
}

// concatenate the data:
process concatenate_gs {
    input:
    val pass
    path file_paths
    val concat_file

    output:
    val concat_file

    script:
    """
    cat ${file_paths} > intermediate-cat.txt
    sort intermediate-cat.txt -r | uniq > ${concat_file}
    rm intermediate-cat.txt
    """
}

// define parameters:
// for subsetting data:

// for converting to h5ad file:
params.h5ad = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad"
params.h5ad_coversion_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/csv-h5ad-conversion.R"

// for computing high variance gene sets:
params.high_fraction_gene_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/high-expression-gene-permutation.R"
gene_numbers = [100, 500, 1000]
params.subsetted_file = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv"
params.gs_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
params.output_dir = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation-95percentile/"

// for computing high fraction:
params.high_variance_gene_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/high-variance-gene-permutation.R"

// for concatenating into a single file:
params.cat_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation-95percentile/concatenated-random-permutation.gs"

// define execution workflow:
workflow {
    high_frac_gs = high_fraction_genes(params.high_fraction_gene_executable, params.subsetted_file, params.gs_file, gene_numbers, params.output_dir)
    high_frac_gs = high_frac_gs.collect()
    high_var_gs = high_variance_genes(high_frac_gs, params.high_variance_gene_executable, params.subsetted_file, params.gs_file, gene_numbers, params.output_dir)
    high_var_gs = high_var_gs.collect()
    // concatenate the generated gs files:
    concat_file_channels = Channel.fromPath("${params.output_dir}seed-*-genes-permuted-weights-*.gs").collect()
    concatenated_gs = concatenate_gs(high_var_gs, concat_file_channels, params.cat_file)
}
