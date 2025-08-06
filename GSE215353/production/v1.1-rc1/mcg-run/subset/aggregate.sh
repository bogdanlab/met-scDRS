#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=5:00:00,h_data=50G
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

current_date=$(date +"%Y-%m-%d")
# aggregate KNN:
Rscript /u/home/l/lixinzhe/project-github/met-scDRS/met-scDRS-method/archive/version-3.0/knn-aggregates.R \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv' \
    --fraction "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-75K-subset-mcg-v1_1_1_rc1.csv" \
    --dr1 'UMAP_1' \
    --dr2 'UMAP_2' \
    --k 5 \
    --threads 5 \
    --output "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/${current_date}-mcg-subset-knn-aggregate.csv"
