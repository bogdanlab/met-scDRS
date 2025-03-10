#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h='!n6454',h_rt=120:00:00,h_data=30G,highp
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
nextflow run /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-1.0/met-scDRS-v1.0-mouse.nf \
    -c /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/met-scDRS-v1.0-run/met-scDRS-v1.0-mcg-74-traits.config

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
