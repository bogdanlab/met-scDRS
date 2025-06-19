# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/GSE215353/revision/v1.1-rc-subset/runs/no_normalization_submission.sh"

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
        "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v1.2/processed-75K-subset-mch.h5ad" \
        "${gs_file}" \
        "mean_var" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var/sampling/' \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var/'

    # treat the cluster nicely:
    sleep 1

done