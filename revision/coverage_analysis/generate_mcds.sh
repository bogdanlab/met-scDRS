#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -j y
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -N makemcds
#$ -l h_rt=23:59:00,h_data=6G
#$ -pe shared 24
# usage() { echo "Usage: bash $0 -a allc_list_file -o output_prefix -b bin_size -h current_directory_path" 1>&2; exit 1; }

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `

if [ "$#" -ne 2 ]; then
    echo "Usage: qsub generate_mcds.sh <allc_table> <output_path>"
    exit 1
fi

allc_table="$1"
output_path="$2"

echo "ALLC table: $allc_table"
echo "Output path: $output_path"
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh

module load anaconda3
#source activate snmcseq
conda activate allcools

# head -n 3 /u/project/jflint/heffel/BICAN3/Tien10k/allc_table.tsv > /u/home/l/lixinzhe/project-geschwind/port/scDRS/tester/allc_table.tsv

mkdir -p "$(dirname "$output_path")"

allcools generate-dataset \
--allc_table "$allc_table" \
--output_path "$output_path" \
--chrom_size_path /u/project/cluo/heffel/BICAN/ref/chromsizes.tsv \
--obs_dim cell \
--cpu 16 \
--chunk_size 50 \
--regions gene_body /u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/feature_out/gene_body.bed \
--quantifiers gene_body count CGN,CHN


# --quantifiers chrom5k hypo-score CGN cutoff=0.9
#--regions chrom5k 5000 \
##--output_path /u/project/cluo/heffel/snm3Cseq_maptest/4brainreg_02192024.mcds \
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo "Done"
