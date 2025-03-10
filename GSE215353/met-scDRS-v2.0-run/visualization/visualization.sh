current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-subset-fraction-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"

# plot for MCG
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-subset-fraction-mcg-cell-type-significance-proportion.png" \
    --plot_type "proportion"

### spatial disection association plot ############################################################
## test if spatial disection pattern can predict risk score in mcg:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 0.1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-mcg"

## test if spatial disection pattern can predict risk score in mch:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 0.1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-mch"

## test if spatial disection pattern can predict risk score in mcg:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-mcg-all"

## test if spatial disection pattern can predict risk score in mch:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/annotation-disease-score-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 1 \
    --num_cutoff 0 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-mch-all"

## Aggregation:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-aggregate-sig-only/"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-aggregate-all/"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --disease "UKB_460K.biochemistry_TotalProtein" \
    --p_cutoff 1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-UKB_460K.biochemistry_TotalProtein-aggregate-all.csv"

### umap plot ##################################################################################### 
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/v2.0-mch-GSE215353-umap-plot/"

# umap plot: 
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/v2.0-mcg-GSE215353-umap-plot/"

### 3D rendering plot #############################################################################
# Define the prefix to remove
prefix="/Volumes/misses/research/Bogdan-lab/GSE215353-aggregate-zscore/GSE215353-aggregate-all/"
for file in "$prefix"*.csv; do
    filename="${file#$prefix}"
    filename="${filename%.csv}"
    python /Users/tardigrade/Desktop/research/Bogdan-lab/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/spatial-rendering.py "$filename"
done