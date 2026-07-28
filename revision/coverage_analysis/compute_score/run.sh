###########################################################################################
######                                    CpG run                                    ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/CpG-single-cell-run/self_extracted/submission.sh"
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CGN_gene_body_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CGN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CGN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CGN_gene_body_full.cov'
    
    # treat the cluster nicely:
    sleep 1

done

###########################################################################################
######                                    non-CpG run                                ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/CpG-single-cell-run/self_extracted/submission.sh"
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CHN_gene_body_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CHN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CHN/full/' \
        "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_full.cov"
    
    # treat the cluster nicely:
    sleep 1

done

###########################################################################################
######                               Create an intermediate                          ######
###########################################################################################
met_scdrs compute_score \
    --h5ad_file '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CHN_gene_body_QC.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --cov_file "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_full.cov" \
    --intermediate "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/intermediate_merged_071726_CHN_gene_body_QC.pkl" \
    --variance_clip 5 \
    --transformation "arcsine" \
    --h5ad_species human \
    --gs-file /u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/KC_75_traits_split.gs9 \
    --gs_species human \
    --out_folder /u/home/l/lixinzhe/project-geschwind/port/scDRS/tester/ \
    --ctrl_match_opt "mean_var_density" \
    --weight_opt "inv_std" \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score False \
    --diagnostic True \
    --diagnostic_dir /u/home/l/lixinzhe/project-geschwind/port/scDRS/tester/ \
    --verbose True
