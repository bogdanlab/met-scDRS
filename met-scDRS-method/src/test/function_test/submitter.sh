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
weight_opt="$4"
output_dir="$5"

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

# output the remaining gs:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/get-remaining-gs.R \
    --scDRS_dir ${output_dir} \
    --gs_file ${gs_file} \
    --output_gs ${tmp_gs}

# compute the score using the remaining gs:
met_scdrs compute_score \
    --h5ad_file ${h5ad_file} \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --h5ad_species mouse \
    --cov_file '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv' \
    --gs-file ${tmp_gs} \
    --gs_species human \
    --out_folder ${output_dir} \
    --ctrl_match_opt ${control_scheme} \
    --weight_opt ${weight_opt} \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic False \
    --diagnostic_dir '/u/scratch/l/lixinzhe/revision_scratch/batch_sampling/' \
    --verbose True

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####