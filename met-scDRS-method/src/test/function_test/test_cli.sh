met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --cov_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv' \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/' \
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score True \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --verbose True
    
met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs' \
    --gs_species human \
    --out_folder '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/' \
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/ \
    --verbose True


###########################################################################################
######                                    parallel                                   ######
###########################################################################################
# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/met-scDRS-method/src/test/function_test/submitter.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}KC_75_traits_split.gs"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad" \
        "${gs_file}" \
        "mean_var" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/mean_var/'

    # treat the cluster nicely:
    sleep 1

done

# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/met-scDRS-method/src/test/function_test/submitter.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}KC_75_traits_split.gs"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad" \
        "${gs_file}" \
        "mean_var_length" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/mean_var_length/'

    # treat the cluster nicely:
    sleep 1

done

###########################################################################################
######                               distribution check                              ######
###########################################################################################
met_scdrs probe_background \
    --score /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/ \
    --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_10K/ \
    --sampling 5 \
    --cell_meta_path /u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-10K-subset-GSE132489-meta.tsv \
    --group_column MajorType \
    --seed 103