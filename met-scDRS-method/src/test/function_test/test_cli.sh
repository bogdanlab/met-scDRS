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

# normalization test
met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --transformation 'library' \
    --h5ad_species mouse \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/mean_var_length/library/' \
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/diagnostic/ \
    --verbose True

###########################################################################################
######                                    Diagnostics                                ######
###########################################################################################
met_scdrs compare_score \
    --score1 "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var/" \
    --score2 "/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/" \
    --plot_path "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K_comparison"

met_scdrs compare_score \
    --score1 "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var/" \
    --score2 "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/" \
    --plot_path "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K_comparison/mean_var_vs_mean_var_length/"

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
    --seed 103 \
    --thread_num 2
    
###########################################################################################
######                                    check save                              ######
###########################################################################################
met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --transformation arcsine \
    --h5ad_species mouse \
    --gs-file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/scratch/l/lixinzhe/tmp_loc/' \
    --intermediate '/u/scratch/l/lixinzhe/tmp_loc/small_test_processed.pkl'\
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score True \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir '/u/scratch/l/lixinzhe/tmp_loc/sampling/' \
    --verbose True