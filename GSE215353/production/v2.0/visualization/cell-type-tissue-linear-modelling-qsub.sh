#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=20G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 3
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
file="$1"
result_dir="$2"

echo "trait file: $file"
echo "output directory: $result_dir"

mkdir ${result_dir}
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/monte-carlo-linear-model.R \
    --score_file ${file} \
    --meta_data "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv" \
    --field1 "X_MajorType" \
    --field2 "tissue" \
    --p_cutoff 0.1 \
    --min_sig_cell 0 \
    --min_cell_num 100 \
    --out "${result_dir}"

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####