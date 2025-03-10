current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-production-fraction-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"

# plot for MCG
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-production-fraction-mcg-cell-type-significance-proportion.png" \
    --plot_type "proportion"

### count heatmap plot ############################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-production-fraction-mch-cell-type-significance-count.png" \
    --plot_type "count"

# plot for MCG
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-production-fraction-mcg-cell-type-significance-count.png" \
    --plot_type "count"

### UMAP visualization ###########################################################################
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/v2.0-mch-GSE215353-production-umap-plot/"

# umap plot: 
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/v2.0-mcg-GSE215353-production-umap-plot/"

### spatial rendering #############################################################################
# first we aggregate the information:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-production-aggregate-all/"

# we also aggregate for significant cells only:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-score-aggregation.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field1 "X_MajorType" \
    --field2 "X_Region" \
    --p_cutoff 0.1 \
    --out "/u/scratch/l/lixinzhe/tmp-file/GSE215353-production-aggregate-sig-only/"

# downloaded  /u/scratch/l/lixinzhe/tmp-file/GSE215353-production-aggregate-all/2025-03-04-PASS_MDD_Howard2019-X_MajorType-X_Regionsummary.csv to /Volumes/misses/research/Bogdan-lab/2024-03-04-region-render/GSE215353-production-MDD-only/
# Define the prefix to remove
prefix="/Volumes/misses/research/Bogdan-lab/2024-03-04-region-render/GSE215353-production-MDD-only/"
for file in "$prefix"*.csv; do
    filename="${file#$prefix}"
    filename="${filename%.csv}"
    python /Volumes/misses/research/Bogdan-lab/2024-03-04-region-render/spatial-rendering.py \
        "$filename" \
        "$prefix" \
        "/Volumes/misses/research/Bogdan-lab/2024-03-04-region-render/GSE215353-production-MDD-only/render/" \
        "true"
done

### mc region modelling ###########################################################################
# compute region's beta:
bash /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/cell-type-tissue-linear-modelling.sh

### VOLCANO PLOT ##################################################################################
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-linear-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv" \
    --field1 "X_MajorType" \
    --field2 "tissue" \
    --p_cutoff 0.1 \
    --num_cutoff 0 \
    --min_cell_num 100 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/"
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/cell-type-tissue-linear-model-distribution.R \
    '/u/home/l/lixinzhe/project-geschwind/plot/2024-04-17-effect-size-X_MajorType-tissue-disease-table.csv'

# plot out the results using all the cells:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/brain-region-linear-model.R \
    --score_dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv" \
    --field1 "X_MajorType" \
    --field2 "tissue" \
    --p_cutoff 1 \
    --num_cutoff 0 \
    --min_cell_num 100 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/"
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/cell-type-tissue-linear-model-distribution.R \
    '/u/home/l/lixinzhe/project-geschwind/plot/2024-04-16-effect-size-X_MajorType-tissue-disease-table.csv'
