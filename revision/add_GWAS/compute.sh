#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=10G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 7
# Email address to notify
#$ -M lxzjason@gmail.com
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

met_scdrs compute_score \
    --h5ad_file '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.pkl' \
    --preprocess False \
    --preprocess_method inverse \
    --variance_clip 5 \
    --transformation arcsine \
    --h5ad_species human \
    --cov_file '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/full-mch-centered-log-library.cov' \
    --gs-file "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-output.gs" \
    --gs_species human \
    --out_folder /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/met_scdrs/out/ \
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score False \
    --diagnostic False \
    --verbose True
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "