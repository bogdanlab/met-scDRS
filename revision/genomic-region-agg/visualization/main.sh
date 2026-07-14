### proportion_heatmap.sh #########################################################################
current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CHN/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "L3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-fraction-mch-promoter-cell-type-significance-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_exon_CHN/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "L3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-fraction-mch-exon-cell-type-significance-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/met-scDRS/revision/genomic-region-agg/visualization/proportion_auto.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_intron_CHN/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv' \
    --field "L3" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revision-GSE215353-50k-fraction-mch-cell-intron-type-significance-proportion.png" \
    --plot_type "proportion" \
    --color_scale_max 0.7

