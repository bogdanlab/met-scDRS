### compare-batch-correction.sh ###################################################################
# compare between batch correction:
met_scdrs compare_score \
    --score1 "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var_length/" \
    --score2 "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges215353_75K_subset/mean_var_length_with_cov/" \
    --plot_path "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/ges215353_75K_subset_comparison_mean_var_len_w_wo_cov"
