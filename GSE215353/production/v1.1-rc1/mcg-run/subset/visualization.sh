current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript  /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_mcg_75K_subset_knn/mean_var_length_with_cov/" \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-v1.1rc-GSE215353-subset-fraction-mcg-cell-type-significance-proportion.png" \
    --plot_type "proportion"
