###########################################################################################
######                                    Visualization                              ######
###########################################################################################
current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CHN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-mch-gene-body-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 1

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CGN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-mcg-gene-body-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 1