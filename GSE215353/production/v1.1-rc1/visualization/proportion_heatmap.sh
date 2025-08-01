### proportion_heatmap.sh #########################################################################
current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/significant-cell-visualization-script.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-production-fraction-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"