met-scdrs compute_score \
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
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True
    
met-scdrs compute_score \
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
    