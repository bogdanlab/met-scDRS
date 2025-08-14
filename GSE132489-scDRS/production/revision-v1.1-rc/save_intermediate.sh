#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=50G
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

input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
gs_file=${input_gs_dir}KC_75_traits_split.gs8
# for each of the gs file submit a job:
echo "read gs file:"
echo "$gs_file"

met_scdrs compute_score \
    --h5ad_file '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.h5ad' \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --transformation arcsine \
    --h5ad_species mouse \
    --cov_file '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.cov' \
    --gs-file ${gs_file} \
    --gs_species human \
    --out_folder '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/' \
    --intermediate '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.pkl'\
    --ctrl_match_opt mean_var_length \
    --weight_opt inv_std \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score True \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/sampling/' \
    --verbose True

 echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####