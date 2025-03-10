#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h='!n6454',h_rt=288:00:00,h_data=25G,highp
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 5
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
# the reason that we are referring to met-scDRS-v1.0 is that we want to call methylation scDRS without using covariate
nextflow run /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/met-scDRS-without-covariates.nf \
    -c /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/production/production-met-scDRS-v2.0/mcg-run/met-scDRS-mcg.config \
    -resume
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####