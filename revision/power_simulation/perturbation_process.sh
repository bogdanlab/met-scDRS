# create h5ad files with different perturbation
executable="/u/home/l/lixinzhe/project-github/met-scDRS/revision/power_simulation/perturbation.R"
gs_file="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
data_matrix='/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.csv'
output_dir="/u/scratch/l/lixinzhe/revision_scratch/simulation/fixed-overlap/"

Rscript ${executable} \
    --data_matrix ${data_matrix} \
    --gs_file ${gs_file} \
    --trait "UKB_460K.body_HEIGHTz" \
    --perturbation_effect_start 0.10 \
    --perturbation_effect_end 0.9 \
    --effect_step 0.1 \
    --gene_number 1000 \
    --cell_number 500 \
    --overlap_start 0.5 \
    --overlap_end 0.5 \
    --overlap_step 0.1 \
    --replication 100 \
    --output_dir ${output_dir}

output_dir="/u/scratch/l/lixinzhe/revision_scratch/simulation/fixed-perturbation/"
Rscript ${executable} \
    --data_matrix ${data_matrix} \
    --gs_file ${gs_file} \
    --trait "UKB_460K.body_HEIGHTz" \
    --perturbation_effect_start 1.005 \
    --perturbation_effect_end 1.005 \
    --effect_step 0.001 \
    --gene_number 1000 \
    --cell_number 500 \
    --overlap_start 0.1 \
    --overlap_end 0.5 \
    --overlap_step 0.1 \
    --replication 100 \
    --output_dir ${output_dir}