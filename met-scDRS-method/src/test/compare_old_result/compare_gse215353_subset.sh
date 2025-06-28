
# recreating v1.0 preprocessed run to make sure that things are comparable:
met_scdrs compute_score \
    --h5ad_file '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/subset_randomized_mch_gene_fraction_copy.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/recreate/gse132489_subset_v1/' \
    --ctrl_match_opt mean_var \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic False \
    --verbose True

# compare to /u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/ results

met_scdrs compare_score \
    --score1_path /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/recreate/gse132489_subset_v1/PASS_Alzheimers_Jansen2019.score.gz \
    --score2_path /u/project/geschwind/lixinzhe/scDRS-output/GSE132489/mch/met-scDRS-v1-run/LXZ-74-traits/PASS_Alzheimers_Jansen2019.score.gz \
    --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/alzheimer_comparison_package_vs_preprint.png