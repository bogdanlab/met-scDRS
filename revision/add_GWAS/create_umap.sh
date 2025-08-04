# create umap of the annotation:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/met_scdrs/out/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse215353-full-mch-length-arcsine-regress/umap/"
