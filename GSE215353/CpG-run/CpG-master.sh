current_date=$(date +"%Y-%m-%d")
# aggregate KNN:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-3.0/knn-aggregates.R \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --fraction '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mcg.csv' \
    --dr1 'UMAP_1' \
    --dr2 'UMAP_2' \
    --k 5 \
    --threads 8 \
    --output "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/${current_date}-mcg-subset-knn-aggregate.csv"

# use the aggregated KNN CpG to data:
python "/u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/csv-to-h5ad.py" \
    "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/2025-01-14-mcg-subset-knn-aggregate.csv" \
    "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mcg-knn5-aggregate.h5ad"

# call scDRS:
input_gs_dir="/u/scratch/l/lixinzhe/tmp-file/job-array/"
out_dir="/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v3/mcg/GSE215353-mcg-knn/"
submission_script="/u/home/l/lixinzhe/project-github/scDDS/application/immune-adult-2022/submitter.sh"
h5ad_file='/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mcg-knn5-aggregate.h5ad'

# split data:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}KC_75_traits_split.gs"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${h5ad_file}" \
        "${out_dir}"

    # treat the cluster nicely:
    sleep 1

done
