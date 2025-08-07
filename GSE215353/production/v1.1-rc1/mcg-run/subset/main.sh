### convert into h5ad from the KNN aggregated results
met_scdrs convert_to_h5ad --csv_path "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/2025-08-05-mcg-subset-knn-aggregate.csv"

### additional preprocessing step: devide all values by number of KNN so it is average instead of sum
python /u/home/l/lixinzhe/project-github/met-scDRS/GSE215353/production/v1.1-rc1/mcg-run/subset/pickle_CpG.py

### Compute the score:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/GSE215353/production/v1.1-rc1/mcg-run/subset/submission.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/2025-08-05-mcg-subset-knn-aggregate.pkl" \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_mcg_75K_subset_knn/mean_var_length_with_cov/sampling/' \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_mcg_75K_subset_knn/mean_var_length_with_cov/' \
        "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/subset75k-mcg-centered-log-library.cov"

    # treat the cluster nicely:
    sleep 1

done
