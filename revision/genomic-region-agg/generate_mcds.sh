#!/bin/bash
#$ -o job$TASK_ID
#$ -cwd
#$ -N makemcds
#$ -l h_rt=2:00:00,h_data=6G
#$ -pe shared 24
# usage() { echo "Usage: bash $0 -a allc_list_file -o output_prefix -b bin_size -h current_directory_path" 1>&2; exit 1; }

cd $wd
echo "CWD: $PWD"

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh

module load anaconda3
#source activate snmcseq
conda activate allcools

# head -n 3 /u/project/jflint/heffel/BICAN3/Tien10k/allc_table.tsv > /u/home/l/lixinzhe/project-geschwind/port/scDRS/tester/allc_table.tsv

allcools generate-dataset \
--allc_table /u/home/l/lixinzhe/project-geschwind/port/scDRS/tester/allc_table.tsv \
--output_path /u/project/geschwind/lixinzhe/data/GSE215353/Tien10k/Tien10k_tester_07092026.mcds \
--chrom_size_path /u/project/cluo/heffel/BICAN/ref/chromsizes.tsv \
--obs_dim cell \
--cpu 16 \
--chunk_size 50 \
--regions promoter /u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/feature_out/promoter_2kb_500bp.bed \
--regions intron /u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/feature_out/intron_segments.bed \
--regions exon /u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/feature_out/exon_merged_segments.bed \
--quantifiers promoter count CGN,CHN \
--quantifiers intron count CGN,CHN \
--quantifiers exon count CGN,CHN \
# --quantifiers chrom5k hypo-score CGN cutoff=0.9
#--regions chrom5k 5000 \
##--output_path /u/project/cluo/heffel/snm3Cseq_maptest/4brainreg_02192024.mcds \
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo "Done"
