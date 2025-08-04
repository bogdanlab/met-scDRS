###########################################################################################
######                                  process sumstat                              ######
###########################################################################################
# process summary statistics:
cd /u/scratch/l/lixinzhe/revision_scratch/mdd_gwas/
gunzip pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz
tail -n +55 pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv > mdd_gwas_without_header.tsv

# formatting:
awk 'BEGIN{ FS = OFS = "\t" } {print $3"\t"$8"\t"$13+$14}' "mdd_gwas_without_header.tsv"  > "MDD-ALL-sumstat_N.tsv"
new_header='SNP\tP\tN'
sed -i "1 s/.*/$new_header/" MDD-ALL-sumstat_N.tsv

# run MAGMA to get scores:
mkdir -p /u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"
${magma_dir}/magma \
    --bfile ${magma_dir}/aux/g1000_eur \
    --pval MDD-ALL-sumstat_N.tsv use='SNP,P' ncol='N' \
    --gene-annot /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1.genes.annot \
    --out /u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd

Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-magma-stats.R" \
    /u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd.genes.out \
    ${magma_dir}/aux/NCBI37.3.gene.loc \
    "/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd-munge-input.txt"

# change the traits header:
new_header="gene\tPGC_MDD_2025"
sed -i "1 s/.*/$new_header/" "/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd-munge-input.txt"

# deploy munge gs for the SCZ magma output
scdrs munge-gs \
    --out-file /u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd-munge-output.gs \
    --zscore-file "/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/ALL/out/step2/mdd-munge-input.txt" \
    --weight zscore \
    --n-max 1000 \
    --n_min 50

