#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h='!n6454',h_rt=48:00:00,h_data=16G, highp
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
cd /u/scratch/l/lixinzhe/process/
gs_file="$1"
h5ad_file="$2"
output_dir="$3"
cov_file="$4"
DIR="/u/scratch/l/lixinzhe/tmp-file/tmp-gs/"
tmp_gs="${DIR}${JOB_ID}_tmp.gs"

# get unfinished jobs:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/get-remaining-gs.R \
    --scDRS_dir ${output_dir} \
    --gs_file ${gs_file} \
    --output_gs ${tmp_gs}

# compute:
scdrs compute-score \
    --h5ad-file ${h5ad_file} \
    --h5ad-species "human" \
    --gs-file ${tmp_gs} \
    --gs-species human \
    --out-folder ${output_dir} \
    --cov-file ${cov_file} \
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