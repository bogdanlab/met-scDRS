// invert the fraction:
process methylation_processing {
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

// run the simulated scDRS
process scDRS_run {
    input:
    path h5ad_path
    path gs_file
    val output_dir
    path cov_file

    output:
    val output_dir

    script:
    """
    scdrs compute-score \
        --h5ad-file ${h5ad_path} \
        --h5ad-species human \
        --gs-file ${gs_file} \
        --gs-species human \
        --out-folder ${output_dir} \
        --cov-file "${cov_file}" \
        --flag-filter-data False \
        --flag-raw-count False \
        --n-ctrl 1000 \
        --weight_opt "inv_std" \
        --flag-return-ctrl-raw-score False \
        --flag-return-ctrl-norm-score True
    """
    
}

// define workfow:
workflow {
    //invert fraction:
    processed_met = methylation_processing(params.processing_executable, params.methylation_csv, params.methylation_processed_output)
    // convert the subset csv file into h5ad object:
    h5ad = convert_h5ad(params.h5ad_coversion_executable, processed_met, params.h5ad_output)
    // compute scDRS:
    scDRS_output_dir = scDRS_run(h5ad, params.gene_set_file, params.output_dir, params.cov_file)
}