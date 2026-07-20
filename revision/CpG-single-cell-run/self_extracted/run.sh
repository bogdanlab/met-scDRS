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
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/'
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CGN_gene_body_QC.h5ad' \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/cov/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/cov/' \
        "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CGN_gene_body_centered_log_rowsum.cov"
    
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
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/'
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CHN_gene_body_QC.h5ad' \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/cov/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/cov/' \
        "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_centered_log_rowsum.cov"
    
    # treat the cluster nicely:
    sleep 1

done