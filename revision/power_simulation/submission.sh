#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=2:00:00,h_data=4G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
PATH=$PATH:$HOME/bin:/u/project/pasaniuc/lixinzhe/software/nextflow
export PATH
module load gcc
module load intel
module load java/jdk-11.0.14
module load anaconda3
conda activate default_r_base

## execute command
# Process the CSV fields as separate parameters
h5ad_file="$1"
gs_file="$2"
control_scheme="$3"
transform_scheme="$4"
weight_opt="$5"
diagnostic_dir="$6"
output_dir="$7"

# NOTE, the last (6th argument must be covariate file, otherwise interpretted incorrectly)
cov_file="$8"

DIR="/u/scratch/l/lixinzhe/tmp-file/tmp-gs/"
tmp_gs="${DIR}${JOB_ID}_tmp.gs"

echo "gs input file: $gs_file"
echo "h5ad input file: $h5ad_file"
echo "output directory: $output_dir"
echo "tmp gene set file: ${tmp_gs}"


# perform a check if the data directory exist, if not create it:
if [ ! -d "$DIR" ]; then
  # If the directory doesn't exist, create it
  mkdir -p "$DIR"
  echo "Directory $DIR created."
else
  echo "Directory $DIR already exists."
fi
mkdir -p $diagnostic_dir

# first we will have to modify the trait name to the right name:
new_trait_name=$(basename "$h5ad_file")
mkdir -p /u/scratch/l/lixinzhe/tmp-file/causal-simulation/
renamed_gs="/u/scratch/l/lixinzhe/tmp-file/causal-simulation/${new_trait_name}.gs"
sed "s/UKB_460K.body_HEIGHTz/${new_trait_name}/" ${gs_file} > ${renamed_gs}

# compute the score using the remaining gs:
if [ -n "$cov_file" ]; then
    met_scdrs compute_score \
        --h5ad_file ${h5ad_file} \
        --preprocess True \
        --preprocess_method inverse \
        --variance_clip 5 \
        --transformation ${transform_scheme} \
        --h5ad_species mouse \
        --cov_file ${cov_file} \
        --gs-file ${renamed_gs} \
        --gs_species human \
        --out_folder ${output_dir} \
        --ctrl_match_opt ${control_scheme} \
        --weight_opt ${weight_opt} \
        --n_ctrl 1000 \
        --flag_return_ctrl_raw_score False \
        --flag_return_ctrl_norm_score False \
        --diagnostic False \
        --diagnostic_dir ${diagnostic_dir} \
        --verbose True
else
    met_scdrs compute_score \
        --h5ad_file ${h5ad_file} \
        --preprocess True \
        --preprocess_method inverse \
        --variance_clip 5 \
        --transformation ${transform_scheme} \
        --h5ad_species mouse \
        --gs-file ${renamed_gs} \
        --gs_species human \
        --out_folder ${output_dir} \
        --ctrl_match_opt ${control_scheme} \
        --weight_opt ${weight_opt} \
        --n_ctrl 1000 \
        --flag_return_ctrl_raw_score False \
        --flag_return_ctrl_norm_score False \
        --diagnostic False \
        --diagnostic_dir ${diagnostic_dir} \
        --verbose True
fi

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####