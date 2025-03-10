submission_script="/u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0-with-cov/mch-run/met-scDRS-mch.sh"
input_gs_dir="/u/scratch/l/lixinzhe/tmp-file/job-array/"
submission_script="/u/home/l/lixinzhe/project-github/scDDS/application/immune-adult-2022/submitter.sh"
h5ad_file='/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mch.h5ad'
out_dir="/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-dev-with-cov/"
cov_file="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-covarite-subset.cov"

# split data:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}GSE215353-mch-with-cov-dev.gs"

for gs_file in ${input_gs_dir}GSE215353-mch-with-cov-dev.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${h5ad_file}" \
        "${out_dir}" \
        "${cov_file}"

    # treat the cluster nicely:
    sleep 1

done
