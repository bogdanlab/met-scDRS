#!/usr/bin/env nextflow
/*
 * pipeline for methylation data processing:
 */

// define subsetting process:
process subset_data {
    input:
    path randomized_data
    val lines
    val subset_file

    output:
    val subset_file

    script:
    """
    head -n ${lines} ${randomized_data} > ${subset_file}
    """
}

// invert the fraction:
process invert_fraction {
    input:
    path executable
    path subset_file
    val inverted_file

    output:
    val inverted_file

    script:
    """
    Rscript ${executable} \
        --data_matrix ${subset_file} \
        --output_file ${inverted_file}
    """
}

// convert the subsetted file to h5ad:
process convert_h5ad {
    input:
    path executable
    path subset_file
    val h5ad_file

    output:
    val h5ad_file

    script:
    """
    Rscript ${executable} \
        --data_matrix ${subset_file} \
        --output_file ${h5ad_file}
    """
}

// define getting the putative gene set from the data:
process random_genes {
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
            --replication 100 \
            --output_dir ${output_dir} > "gene_num_${gene_number}.logfile"

    """
}

// define getting the high fraction gene set:
process high_fraction_genes {
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
            --quantile 0.75 \
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
            --quantile 0.75 \
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

// run the simulated scDRS
process scDRS_run {
    input:
    path h5ad_path
    path gs_file
    val output_dir

    output:
    val output_dir

    script:
    """
    scdrs compute-score \
        --h5ad-file ${h5ad_path} \
        --h5ad-species mouse \
        --gs-file ${gs_file} \
        --gs-species mouse \
        --out-folder ${output_dir} \
        --flag-filter-data False \
        --flag-raw-count False \
        --n-ctrl 1000 \
        --flag-return-ctrl-raw-score False \
        --flag-return-ctrl-norm-score True
    """
    
}


// define parameters:
// for subsetting data:
params.lines = 10001
params.full_fraction_file = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.csv"
params.subset_fraction_output = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.csv"

// for inverting the fraction:
params.inverted_fraction = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv"
params.invert_fraction_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/inverse-fraction.R"

// for converting to h5ad file:
params.h5ad = "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad"
params.h5ad_coversion_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/csv-h5ad-conversion.R"

// for computing random gene sets:
gene_numbers = [100, 500, 1000]
params.random_gene_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/random-genes-permutation.R"
params.gs_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
params.output_dir = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/"

// for computing high fraction gene sets:
params.high_fraction_gene_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/high-expression-gene-permutation.R"

// for computing high fraction:
params.high_variance_gene_executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/high-variance-gene-permutation.R"

// for concatenating into a single file:
params.cat_file = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/concatenated-random-permutation.gs"

// for computing scDRS:
params.scdrs_output_dir = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/random-simulation-output/"

// define execution workflow:
workflow {
    // subset permuted data:
    subsetted_file = subset_data(params.full_fraction_file, params.lines, params.subset_fraction_output)
    //invert fraction:
    inverted_fraction = invert_fraction(params.invert_fraction_executable, subsetted_file, params.inverted_fraction)
    // convert the subset csv file into h5ad object:
    h5ad = convert_h5ad(params.h5ad_coversion_executable, inverted_fraction, params.h5ad)
    // get permuted gene set file for running null simulation
    random_gs = random_genes(h5ad, params.random_gene_executable, inverted_fraction, params.gs_file, gene_numbers, params.output_dir)
    random_gs = random_gs.collect()
    high_frac_gs = high_fraction_genes(random_gs, params.high_fraction_gene_executable, inverted_fraction, params.gs_file, gene_numbers, params.output_dir)
    high_frac_gs = high_frac_gs.collect()
    high_var_gs = high_variance_genes(high_frac_gs, params.high_variance_gene_executable, inverted_fraction, params.gs_file, gene_numbers, params.output_dir)
    high_var_gs = high_var_gs.collect()
    // concatenate the generated gs files:
    concat_file_channels = Channel.fromPath("${params.output_dir}seed-*-genes-permuted-weights-*.gs").collect()
    concatenated_gs = concatenate_gs(high_var_gs, concat_file_channels, params.cat_file)
    // compute scDRS:
    // scDRS_output_dir = scDRS_run(h5ad, concatenated_gs, params.scdrs_output_dir)
}
