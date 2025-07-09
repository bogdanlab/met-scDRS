### scz_gs_score.sh ###############################################################################
# specify common file paths:
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"
mkdir -p /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1

# specify input file:
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/aux/g1000_eur.bim \
    --gene-loc ${magma_dir}/aux/NCBI37.3.gene.loc \
    --out /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1

###########################################################################################
######                                  process sumstat                              ######
###########################################################################################
# process summary statistics:
tail -n +74 PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv > SCZ-EUR-sumstat.tsv
awk 'BEGIN{ FS = OFS = "\t" } {print $2"\t"$11"\t"$12+$13}' "SCZ-EUR-sumstat.tsv" > "SCZ-EUR-sumstat_N.tsv"
new_header='SNP\tP\tN'
sed -i "1 s/.*/$new_header/" SCZ-EUR-sumstat_N.tsv

# run MAGMA to get scores:
mkdir -p /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2
${magma_dir}/magma \
    --bfile ${magma_dir}/aux/g1000_eur \
    --pval /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/SCZ-EUR-sumstat_N.tsv use='SNP,P' ncol='N' \
    --gene-annot /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1.genes.annot \
    --out /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz

Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-magma-stats.R" \
    "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz.genes.out" \
    ${magma_dir}/aux/NCBI37.3.gene.loc \
    "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz-munge-input.txt"

# change the traits header:
new_header="gene\tTrubetskoy_pardinas_2022"
sed -i "1 s/.*/$new_header/" "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz-munge-input.txt"

# deploy munge gs for the SCZ magma output
scdrs munge-gs \
    --out-file /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz-munge-output.gs \
    --zscore-file "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step2/scz-munge-input.txt" \
    --weight zscore \
    --n-max 1000 \
    --n_min 50
