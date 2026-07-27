
###########################################################################################
######                                    non-CpG run                                ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/met-scDRS/revision/gene_validation/submission.sh"
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

for pickled_batch in /u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/preprocessed_pkl_gene_batch/*.pkl; do
    # get the output batch
    batch_name=$(basename "$pickled_batch" _data.pkl)
    out_dir="/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/batch/${batch_name}/"
    mkdir -p "$out_dir"
    mkdir -p "${out_dir}/sampling/"
    
    for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
        # for each of the gs file submit a job:
        echo "read gs file:"
        echo "$gs_file"
        
        # get expected output file:
        trait_name=$(awk 'NR==2 {print $1}' "$gs_file")
        expected_score_file="${out_dir}/${trait_name}.score.gz"

        if [[ -s "$expected_score_file" ]]; then
            echo "skip existing output:"
            echo "$expected_score_file"
        else
            echo "submit gs file:"
            echo "$gs_file"
            echo "expected output:"
            echo "$expected_score_file"
            
            # compute scDRS:
            qsub ${submission_script} \
                "${pickled_batch}" \
                "${gs_file}" \
                "mean_var_density" \
                "arcsine" \
                "inv_std" \
                "${out_dir}/sampling/" \
                "${out_dir}" \
                "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_full.cov"
            # treat the cluster nicely:
            sleep 1
        fi
    done
done