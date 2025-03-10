## execute command
## first we will compute things for the random simulation:
# define input parameters:
input_gs_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/"
submission_script="/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-simulation/simulation-scdrs-parallel-submission.sh"
mch_h5ad="/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad"

# from the data, grab out the gs files and submit job:
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/OCT22-2024-ivd/RANDOM/"
for gs_file in ${input_gs_dir}seed-*-random-genes-permuted-weights-*.gs; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${mch_h5ad}" \
        "${output_dir}"

    # treat the cluster nicely:
    sleep 1

done

## next we will compute things for the high fraction simulation:
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/OCT22-2024-ivd/TOP_FRACTION/"
for gs_file in ${input_gs_dir}seed-*-high-fraction-genes-permuted-weights-*.gs; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${mch_h5ad}" \
        "${output_dir}"

    # treat the cluster nicely:
    sleep 1

done

## next we will compute things for the high variance simulation:
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/OCT22-2024-ivd/TOP_FRACTION/TOP_VARIANCE/"
for gs_file in ${input_gs_dir}seed-*-high-variance-genes-permuted-weights-*.gs; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${mch_h5ad}" \
        "${output_dir}"

    # treat the cluster nicely:
    sleep 1

done
