Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-production-aggregate-all/"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-production-aggregate-sig/"
