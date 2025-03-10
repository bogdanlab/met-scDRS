# define input parameters:
input_gs_file="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/UKB_460K.body_HEIGHTz_only.gs"
submission_script="/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-causal-simulation/simulation-scdrs-parallel-submission.sh"

# fixed overlap
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/result/rare-cell-type/"
mch_h5ad_dir="/u/scratch/l/lixinzhe/tmp-file/causal-simulation/rare-cell-type/"

# first submit the data for the VLMC cell type:
# from the data, grab out the h5ad files and submit job for fixed overlap:
for mch_h5ad in ${mch_h5ad_dir}seed-*-cell-type-*-down-sample-*-proportion-0.25-effect-1.005-overlap-0.25-causal-simulation.h5ad; do
    # for each of the gs file submit a job:
    echo "read mch file:"
    echo "$mch_h5ad"

    # compute scDRS:
    qsub ${submission_script} \
        "${input_gs_file}" \
        "${mch_h5ad}" \
        "${output_dir}"

    # treat the cluster nicely:
    sleep 1

done

# for mch_h5ad in ${mch_h5ad_dir}seed-*-cell-type-ASC-down-sample-*-proportion-0.25-effect-1.005-overlap-0.25-causal-simulation.h5ad; do
#     # for each of the gs file submit a job:
#     echo "read mch file:"
#     echo "$mch_h5ad"

#     # compute scDRS:
#     qsub ${submission_script} \
#         "${input_gs_file}" \
#         "${mch_h5ad}" \
#         "${output_dir}"

#     # treat the cluster nicely:
#     sleep 1

# done
