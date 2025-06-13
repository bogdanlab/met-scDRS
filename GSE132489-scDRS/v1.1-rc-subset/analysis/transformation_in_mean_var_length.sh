# compare betweeen the 30K subset between mean var length and mean var:


# difference between mean var length without vs without normalization
# difference using arcsine ~0.9 spearman correlation
met_scdrs compare_score  \
    --score1 /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/ \
    --score2  /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_arcsine/ \
    --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/arcsine/

# library normalization ~0.99
met_scdrs compare_score  \
    --score1 /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/ \
    --score2  /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_library/ \
    --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/library/

# logit normalization ~ 0.5 - 0.7
met_scdrs compare_score  \
    --score1 /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/ \
    --score2 /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_logit/ \
    --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/logit/
