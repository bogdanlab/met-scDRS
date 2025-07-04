import subprocess

# invoke subprocess:
for gene_num in [100, 500, 1000]:
    # high fraction
    cmd = f"""\
    Rscript /u/home/l/lixinzhe/project-github/met-scDRS/scDRS-methylation-simulation/high-expression-gene-permutation.R \
        --data_matrix '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/normalized-preprocessed-simulation-subset-GSE132489-mch.csv' \
        --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs" \
        --gene_number "{gene_num}" \
        --quantile "0.75" \
        --replication "100" \
        --output_dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-75percentile/"
    """
    print(f'getting gene set for null simulation with number of genes {gene_num}')
    subprocess.run(cmd, shell=True)
    
    # high variance:
    cmd = f"""\
    Rscript  /u/home/l/lixinzhe/project-github/met-scDRS/scDRS-methylation-simulation/high-variance-gene-permutation.R \
        --data_matrix '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/normalized-preprocessed-simulation-subset-GSE132489-mch.csv' \
        --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs" \
        --gene_number "{gene_num}" \
        --quantile "0.75" \
        --replication "100" \
        --output_dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-75percentile/high_variance/"
    """
    subprocess.run(cmd, shell=True)
