# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/power_simulation/submission.sh"

# define input parameters:
input_gs_file="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/UKB_460K.body_HEIGHTz_only.gs"

###########################################################################################
######                                    Fixed overlap                              ######
###########################################################################################
# fixed overlap
output_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/power_simulation/fixed-overlap/informed/"
mch_h5ad_dir="/u/scratch/l/lixinzhe/revision_scratch/simulation/fixed-overlap/informed/"
mkdir -p $output_dir

# from the data, grab out the h5ad files and submit job for fixed overlap:
effect_sizes=("0.6" "0.7" "0.8" "0.9")
for size in "${effect_sizes[@]}"
do
    for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-${size}-overlap-0.5-causal-simulation.h5ad; do
        # for each of the gs file submit a job:
        echo "read mch file:"
        echo "$mch_h5ad"
        basic_name=$(basename "$mch_h5ad")
        output_file="$output_dir$basic_name.score.gz"

        # add an if statement to check if the output already exist:
        if [ -e "$output_file" ]; then
            echo "The file $output_file exists."
        else
        # compute scDRS:
        qsub ${submission_script} \
            "${mch_h5ad}" \
            "${input_gs_file}" \
            "mean_var_length" \
            "arcsine" \
            "inv_std" \
            ${output_dir} \
            ${output_dir}

            # treat the cluster nicely:
            sleep 1
        fi
    done
done

###########################################################################################
######                                    Fixed effect                               ######
###########################################################################################
# fixed effect:
# output_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/power_simulation/fixed-perturbation/"
# mch_h5ad_dir="/u/scratch/l/lixinzhe/revision_scratch/simulation/fixed-perturbation/"
# mkdir -p $output_dir

# # from the data, grab out the h5ad files and submit job for fixed overlap:
# for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-1.005-overlap-*-causal-simulation.h5ad; do
#     # for each of the gs file submit a job:
#     echo "read mch file:"
#     echo "$mch_h5ad"
#     basic_name=$(basename "$mch_h5ad")
#     output_file="$output_dir$basic_name.score.gz"

#     # add an if statement to check if the output already exist:
#     if [ -e "$output_file" ]; then
#         echo "The file $output_file exists."
#     else
#         # compute scDRS:
#         qsub ${submission_script} \
#             "${mch_h5ad}" \
#             "${input_gs_file}" \
#             "mean_var_length" \
#             "arcsine" \
#             "inv_std" \
#             ${output_dir} \
#             ${output_dir}


#         # treat the cluster nicely:
#         sleep 1
#     fi
# done