### umap.sh #######################################################################################
# purpose: visualize umap for methylation on MDD

### PROCESS #######################################################################################
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse215353-full-mch-length-arcsine-regress/umap/"
