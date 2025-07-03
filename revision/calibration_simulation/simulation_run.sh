# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/calibration_simulation/simulation_submission.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-95percentile/"

for gs_file in ${input_gs_dir}*gs; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        "/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad" \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-95percentile/result/' \
        '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-95percentile/result/'

    # treat the cluster nicely:
    sleep 1

done