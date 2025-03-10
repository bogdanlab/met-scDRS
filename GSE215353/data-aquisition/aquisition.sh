### aquisition.sh #################################################################################
# purpose: serves as a master script that will get the data into a csv format for met-scDRS

### AQUISITION ####################################################################################
# extract data out of the rds object for gene level methylation 

## MCH:
# for non neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/e39765a9-3d2c-4296-9f78-26fc38d4a4a0.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/"

# for inhibitory neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/940105e9-6d36-4a70-9712-de3a808c4a09.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/"

# for excitatory neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/606c03ae-c591-457f-849a-6d29ee5390e5.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/"

## MCG:
# for non neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mcg/6bc0a805-a1b2-4c0d-8891-cb5c8804a927.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/"

# for inhibitory neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mcg/74eacfbe-bec4-4348-b9b4-be17630baab8.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/"

# for excitatory neurons:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/data-aquisition/data-extraction.R \
    --rds "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mcg/93ea50ca-85ac-48fb-bc1f-634c247d6d64.rds" \
    --percentage 0.2 \
    --out "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/"

### CONCATENATE ###################################################################################
### MCH:
# grab out the file 1 2 and 3:
file1="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/e39765a9-3d2c-4296-9f78-26fc38d4a4a0-fraction.csv"
file2="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/940105e9-6d36-4a70-9712-de3a808c4a09-fraction.csv"
file3="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/606c03ae-c591-457f-849a-6d29ee5390e5-fraction.csv"

# if the genes are the same order, concatenate the two files:
# Read the first line of each file
first_line_file1=$(head -n 1 "$file1")
first_line_file2=$(head -n 1 "$file2")
first_line_file3=$(head -n 1 "$file3")

# Compare the first lines
if [ "$first_line_file1" = "$first_line_file2" ] && [ "$first_line_file1" = "$first_line_file3" ]; then
    echo "The first lines of $file1, $file2, and $file3 are the same."
    head -n1 $file1 > /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/GSE215353-concat.csv
    tail -n +2 $file1 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/GSE215353-concat.csv
    tail -n +2 $file2 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/GSE215353-concat.csv
    tail -n +2 $file3 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/GSE215353-concat.csv
else
    echo "The first lines of $file1, $file2, and $file3 are different."
fi

### MCG:
file1="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/6bc0a805-a1b2-4c0d-8891-cb5c8804a927-fraction.csv"
file2="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/74eacfbe-bec4-4348-b9b4-be17630baab8-fraction.csv"
file3="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/93ea50ca-85ac-48fb-bc1f-634c247d6d64-fraction.csv"

# if the genes are the same order, concatenate the two files:
# Read the first line of each file
first_line_file1=$(head -n 1 "$file1")
first_line_file2=$(head -n 1 "$file2")
first_line_file3=$(head -n 1 "$file3")

# Compare the first lines
if [ "$first_line_file1" = "$first_line_file2" ] && [ "$first_line_file1" = "$first_line_file3" ]; then
    echo "The first lines of $file1, $file2, and $file3 are the same."
    head -n1 $file1 > /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/GSE215353-concat.csv
    tail -n +2 $file1 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/GSE215353-concat.csv
    tail -n +2 $file2 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/GSE215353-concat.csv
    tail -n +2 $file3 >> /u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/GSE215353-concat.csv
else
    echo "The first lines of $file1, $file2, and $file3 are different."
fi
