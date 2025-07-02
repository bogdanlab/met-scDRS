### submission.sh #################################################################################
# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/GSE215353/production/v1.1-rc1/mch-run/submission.sh"

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
        '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.pkl' \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/' \
        '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/full-mch-centered-log-library.cov'

    # treat the cluster nicely:
    sleep 1

done
