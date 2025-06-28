#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=3:00:00,h_data=7G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-4

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

case $SGE_TASK_ID in
1)
    # No transformation:
    met_scdrs probe_background \
        --score /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/ \
        --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/untransformed/null_distribution/ \
        --sampling 5 \
        --cell_meta_path /u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv \
        --group_column MajorType \
        --seed 103
    ;;
2)
    # arcsine:
    met_scdrs probe_background \
        --score /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_arcsine/ \
        --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/arcsine/null_distribution/ \
        --sampling 5 \
        --cell_meta_path /u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv \
        --group_column MajorType \
        --seed 103
    ;;
3)
    # logit:
    met_scdrs probe_background \
        --score /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_logit/ \
        --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/logit/null_distribution/ \
        --sampling 5 \
        --cell_meta_path /u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv \
        --group_column MajorType \
        --seed 103
    ;;

4)
    # library size
    met_scdrs probe_background \
        --score /u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length_library/ \
        --plot_path /u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/library/null_distribution/ \
        --sampling 5 \
        --cell_meta_path /u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv \
        --group_column MajorType \
        --seed 103
esac
