# first grab out the causal gene set file which is the height trait from the 74 trait gs file:
gs_path="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
head -n1 ${gs_path} > "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/UKB_460K.body_HEIGHTz_only.gs"
grep -E 'UKB_460K.body_HEIGHTz' ${gs_path} >> "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/UKB_460K.body_HEIGHTz_only.gs"

## first we will compute things scDRS:
# define input parameters:
input_gs_file="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/UKB_460K.body_HEIGHTz_only.gs"
submission_script="/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-methylation-causal-simulation/simulation-scdrs-parallel-submission.sh"

# fixed effect:
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/result/fixed-effect/"
mch_h5ad_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-effect/"
# from the data, grab out the h5ad files and submit job for fixed overlap:
for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-1.005-overlap-*-causal-simulation.h5ad; do
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
            "${input_gs_file}" \
            "${mch_h5ad}" \
            "${output_dir}"

        # treat the cluster nicely:
        sleep 1
    fi
done

# fixed overlap
output_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/result/fixed-overlap/"
mch_h5ad_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-overlap/"

# from the data, grab out the h5ad files and submit job for fixed overlap:
effect_sizes=("1.001" "1.002" "1.003" "1.004")
for size in "${effect_sizes[@]}"
do
    for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-${size}-overlap-0.25-causal-simulation.h5ad; do
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
                "${input_gs_file}" \
                "${mch_h5ad}" \
                "${output_dir}"

            # treat the cluster nicely:
            sleep 1
        fi
    done
done

effect_sizes=("1.005" "1.006" "1.007" "1.008")
for size in "${effect_sizes[@]}"
do
    for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-${size}-overlap-0.25-causal-simulation.h5ad; do
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
                "${input_gs_file}" \
                "${mch_h5ad}" \
                "${output_dir}"

            # treat the cluster nicely:
            sleep 1
        fi
    done
done

effect_sizes=("1.009" "1.01" "1.02" "1.03")
for size in "${effect_sizes[@]}"
do
    for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-${size}-overlap-0.25-causal-simulation.h5ad; do
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
                "${input_gs_file}" \
                "${mch_h5ad}" \
                "${output_dir}"

            # treat the cluster nicely:
            sleep 1
        fi
    done
done

effect_sizes=("1.04" "1.05")
for size in "${effect_sizes[@]}"
do
    for mch_h5ad in ${mch_h5ad_dir}seed-*-effect-${size}-overlap-0.25-causal-simulation.h5ad; do
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
                "${input_gs_file}" \
                "${mch_h5ad}" \
                "${output_dir}"

            # treat the cluster nicely:
            sleep 1
        fi
    done
done