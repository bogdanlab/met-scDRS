// get jobs that are not finished yet:
process get_unfinished_jobs {
    input:
    path executable
    path scdrs_folder
    path gs_file
    val remaining_gs_file

    output:
    val remaining_gs_file

    script:
    """
    Rscript ${executable} \
        --scDRS_dir ${scdrs_folder} \
        --gs_file ${gs_file} \
        --output_gs ${remaining_gs_file}
    """
}

// submit jobs again:
process scDRS_run {
    input:
    path h5ad_path
    path remaining_gs_file
    val output_dir
    val species

    output:
    val output_dir

    script:
    """
    scdrs compute-score \
        --h5ad-file ${h5ad_path} \
        --h5ad-species ${species} \
        --gs-file ${remaining_gs_file} \
        --gs-species human \
        --out-folder ${output_dir} \
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
    //get unfinished jobs:
    unfinished = get_unfinished_jobs(params.get_unfinished_executable, params.scdrs_folder, params.gs_file, params.remaining_gs_file)
    // submit those jobs that are not finished:
    scDRS_output_dir = scDRS_run(params.h5ad, unfinished, params.scdrs_folder, params.species)
}