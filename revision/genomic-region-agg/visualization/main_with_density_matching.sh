### proportion_heatmap.sh #########################################################################
current_date=$(date +"%Y-%m-%d")

###########################################################################################
######                    Visualize single cell non-CpG with QC                      ######
###########################################################################################
Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CHN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mch-exon-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CHN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mch-promoter-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CHN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mch-intron-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5

###########################################################################################
######                      Visualize single cell CpG with QC                        ######
###########################################################################################
Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_exon_CGN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mcg-exon-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_promoter_CGN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mcg-promoter-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_intron_CGN/full/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "newL3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-qc-regression-mcg-intron-cell-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.5
