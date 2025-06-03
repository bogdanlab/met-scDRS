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
        --flag-filter-data True \
        --flag-raw-count True \
        --n-ctrl 1000 \
        --flag-return-ctrl-raw-score False \
        --flag-return-ctrl-norm-score True
    """
    
}

// define workfow:
workflow {
    // compute scDRS:
    scDRS_output_dir = scDRS_run(params.h5ad_file, params.gene_set_file, params.output_dir, params.cov_file)
}
