# Specify the folder path
folder="/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full"

# Loop over the files in the folder
for file in "$folder"/*.full_score.gz
do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        echo "Processing file: $file"
        # Get the base name of the file without extension
        base_name=$(basename "$file" .full_score.gz)
        # Extract the part before the first dot
        name_part="${base_name%.full_score.gz}"

        # check if the destination folder exist:
        result_dir="/u/scratch/l/lixinzhe/tmp-file/region-cell-type-analysis-sig-only/${name_part}/"
        if [ -d $result_dir ]; then
            echo "Result already computed, closing"
        else
            echo "Processing $name_part"
            qsub /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/cell-type-tissue-linear-modelling-qsub.sh \
                "$file" \
                "$result_dir"
        fi
    fi
done

# visualization:
# Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/production/v2.0/visualization/cell-type-tissue-mc-model-distribution.R \
#     "/u/scratch/l/lixinzhe/tmp-file/mc_test/PASS_MDD_Howard2019/"
