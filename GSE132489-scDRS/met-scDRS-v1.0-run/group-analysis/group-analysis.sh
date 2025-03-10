### group-analysis.sh #############################################################################
# purpose: perform group analysis for GSE 132489:

## heterogeneity test:
# first add the cell type annotation into the h5ad file:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/add-annotation-h5ad.R \
    --fraction "/u/project/geschwind/lixinzhe/data/processed-met-GSE132489-mch.csv" \
    --meta_data "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv" \
    --field "MajorType" \
    --out "/u/project/geschwind/lixinzhe/data/GSE132489/processed-met-GSE132489-mch-with-major-cell-type.h5ad"

# perform downstream heterogeneity analysis for mch:
scdrs perform-downstream \
    --h5ad-file "/u/project/geschwind/lixinzhe/data/GSE132489/processed-met-GSE132489-mch-with-major-cell-type.h5ad" \
    --score-file "/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/@.full_score.gz" \
    --out-folder /u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/heterogeneity-analysis/ \
    --group-analysis annotation \
    --flag-filter-data False \
    --flag-raw-count False

# aggregate the result into a table:
current_date=$(date +"%Y-%m-%d")
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-1.2/visualization/aggregate-group-analysis-table.R \
    --result_dir '/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/heterogeneity-analysis/' \
    --field 'assoc_mcp' \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v1.0-GSE132489-mch-cell-type-group-analysis-assoc-mcp-fdr.csv"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-1.2/visualization/aggregate-group-analysis-table.R \
    --result_dir '/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/heterogeneity-analysis/' \
    --field 'hetero_mcp' \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v1.2-GSE132489-mch-cell-type-group-analysis-hetero_mcp-fdr.csv"

## test if spatial disection pattern can predict risk score:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/" \
    --meta_data "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv" \
    --field1 "RegionName" \
    --field2 "empty" \
    --p_cutoff 0.1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE132489"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/" \
    --meta_data "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv" \
    --field1 "MajorType" \
    --field2 "RegionName" \
    --p_cutoff 0.1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE132489"
