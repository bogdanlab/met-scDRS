# prepare data for the operation

### DATA:
met_scdrs convert_to_h5ad "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.csv"

### get covariate:
python /u/home/l/lixinzhe/project-github/met-scDRS/GSE132489-scDRS/production/revision-v1.1-rc/rowSum_mch.py

### generate intermediate file:
qsub /u/home/l/lixinzhe/project-github/met-scDRS/GSE132489-scDRS/production/revision-v1.1-rc/save_intermediate.sh

### visualize the MDD UMAP:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/' \
    --meta_data '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.csv' \
    --xaxis "L1UMAP_0" \
    --yaxis "L1UMAP_1" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489-full-mch-length-arcsine-regress/umap/"

### for all traits:
# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/GSE132489-scDRS/production/revision-v1.1-rc/submission.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.pkl' \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/sampling/' \
        '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/' \
        '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.cov'

    # treat the cluster nicely:
    sleep 1

done

current_date=$(date +"%Y-%m-%d")

### proportion heatmap plot #######################################################################
# plot for MCH
Rscript /u/home/l/lixinzhe/project-github/met-scDRS/GSE132489-scDRS/v1.1-rc-subset/analysis/visualization/proportion_heatmap_pipeline.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/" \
    --meta_data '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.csv' \
    --modality "full_mean_var_length_arcsine_regress" \
    --p_cutoff 0.1
