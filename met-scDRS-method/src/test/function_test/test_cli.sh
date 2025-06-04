met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --cov_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv' \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/scratch/l/lixinzhe/revision_scratch/v1.0.0-rc1' \
    --ctrl_match_opt mean_var \
    --weight_opt vs \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --verbose True
    
met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs' \
    --gs_species human \
    --out_folder '/u/scratch/l/lixinzhe/revision_scratch/v1.0.0-rc1' \
    --ctrl_match_opt mean_var \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True

###########################################################################################
######                                 all gene set run                              ######
###########################################################################################
# used to get the background for all genes set in mouse:
met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --cov_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv' \
    --gs_file '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs' \
    --gs_species human \
    --out_folder '/u/scratch/l/lixinzhe/revision_scratch/v1.0.0-rc1' \
    --ctrl_match_opt mean_var \
    --weight_opt vs \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir '/u/scratch/l/lixinzhe/revision_scratch/batch_sampling/' \
    --verbose True