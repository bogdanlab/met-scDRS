###########################################################################################
######                                    CpG run                                    ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/submission.sh"
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CGN_exon_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CGN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CGN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CGN_exon_full.cov'
    
    # treat the cluster nicely:
    sleep 1

    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CGN_intron_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CGN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CGN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CGN_intron_full.cov'
    
    # treat the cluster nicely:
    sleep 1
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CGN_promoter_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CGN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CGN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CGN_promoter_full.cov'
    
    # treat the cluster nicely:
    sleep 1
done

###########################################################################################
######                                    Non-CpG RUn                                ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/submission.sh"
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CHN_exon_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CHN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CHN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CHN_exon_full.cov'
    
    # treat the cluster nicely:
    sleep 1

    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CHN_intron_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CHN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CHN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CHN_intron_full.cov'
    
    # treat the cluster nicely:
    sleep 1
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/merged_073126_CHN_promoter_QC.h5ad' \
        "${gs_file}" \
        "mean_var_density" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CHN/full/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CHN/full/' \
        '/u/project/geschwind/lixinzhe/data/GSE215353/feature_agg_with_coverage/merged_QCed/CHN_promoter_full.cov'
    
    # treat the cluster nicely:
    sleep 1
done

# ### submission.sh #################################################################################
# # call scDRS:

# # we can probably try to do promoter with the inverse at single cell resolution.
# # since promoter CpG is repressing
# # we can check aggregation later

# # intron and exon we probably have to try both ways
# submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/non-CpG_submission.sh"

# # split data:
# input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
# Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
#     --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
#     --output_gs "${input_gs_dir}KC_75_traits_split.gs"

# for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
#     # for each of the gs file submit a job:
#     echo "read gs file:"
#     echo "$gs_file"
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CGN_exon_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_exon_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_exon_CGN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CGN_promoter_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CGN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CGN_intron_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_intron_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_intron_CGN/'
    
#     # treat the cluster nicely:
#     sleep 1

# done

# ###########################################################################################
# ######                                    QCed run                                   ######
# ###########################################################################################
# submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/non-CpG_submission.sh"
# input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

# for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
#     # for each of the gs file submit a job:
#     echo "read gs file:"
#     echo "$gs_file"
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CGN_exon_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CGN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CGN_promoter_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CGN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CGN_intron_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CGN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CGN/'
    
#     # treat the cluster nicely:
#     sleep 1

# done

# ### submission.sh #################################################################################
# # call scDRS:
# submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/non-CpG_submission.sh"

# # split data:
# input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
# Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
#     --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
#     --output_gs "${input_gs_dir}KC_75_traits_split.gs"

# for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
#     # for each of the gs file submit a job:
#     echo "read gs file:"
#     echo "$gs_file"
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CHN_exon_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_exon_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_exon_CHN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CHN_promoter_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CHN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged/merged_07132026_CHN_intron_raw.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_intron_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_intron_CHN/'
    
#     # treat the cluster nicely:
#     sleep 1

# done

# ###########################################################################################
# ######                                    QCed run                                   ######
# ###########################################################################################
# submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/met-scDRS-run/non-CpG_submission.sh"
# input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

# for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
#     # for each of the gs file submit a job:
#     echo "read gs file:"
#     echo "$gs_file"
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CHN_exon_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CHN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CHN_promoter_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CHN/'
    
#     # compute scDRS:
#     qsub ${submission_script} \
#         '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed/merged_071526_CHN_intron_QC.h5ad' \
#         "${gs_file}" \
#         "mean_var_length" \
#         "arcsine" \
#         "inv_std" \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CHN/sampling/' \
#         '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CHN/'
    
#     # treat the cluster nicely:
#     sleep 1

# done
