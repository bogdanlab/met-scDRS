#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=00:40:00,h_data=40G
## Modify the parallel environment
## and the number of cores as needed:
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
gs_file="$1"
h5ad_file="$2"
output_dir="$3"

echo "gs input file: $gs_file"
echo "h5ad input file: $h5ad_file"
echo "output directory: $output_dir"

# first we will have to modify the trait name to the right name:
new_trait_name=$(basename "$h5ad_file")
renamed_gs="/u/scratch/l/lixinzhe/tmp-file/causal-simulation/${new_trait_name}.gs"
sed "s/UKB_460K.body_HEIGHTz/${new_trait_name}/" ${gs_file} > ${renamed_gs}

# then we can compute the scDRS score:
scdrs compute-score \
    --h5ad-file "${h5ad_file}" \
    --h5ad-species mouse \
    --gs-file "${renamed_gs}" \
    --gs-species human \
    --out-folder ${output_dir} \
    --flag-filter-data False \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --weight_opt "inv_std" \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####