### proportion-visualization.sh ###################################################################
# purpose: visualize the proportion change among cell types:

current_date=$(date +"%Y-%m-%d")

# mean_var:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var/' \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-75K-subset-mean-var-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"
    
# mean_var_length:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var_length/' \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-75K-subset-mean-var-length-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"

# mean_var_length_arcsine:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var_length_arcsine/' \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-75K-subset-mean-var-length-arcsine-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"

# mean_var_length_library:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_75K_subset/mean_var_length_library/' \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-75K-subset-mean-var-length-library-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"

# mean_var_length_library:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-script.R \
    --dir '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_75K_subset/mean_var_length_logit/' \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv" \
    --field "X_MajorType" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-GSE215353-75K-subset-mean-var-length-logit-mch-cell-type-significance-proportion.png" \
    --plot_type "proportion"