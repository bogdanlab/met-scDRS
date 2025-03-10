### visualization-master.sh #######################################################################
# purpose: master shell script that documents all visualization of the CpG subset run
current_date=$(date +"%Y-%m-%d")

# visualize the proportional heatmap:
# plot for MCG
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/CpG-run/visualization/proportion-visualization.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE215353-subset/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v2.0-GSE215353-subset-fraction-mcg-cell-type-significance-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/CpG-run/visualization/proportion-visualization.R \
    --dir "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v3/mcg/GSE215353-mcg-knn/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v3-GSE215353-subset-knn-mcg-cell-type-significance-proportion.png" \
    --plot_type "proportion"
