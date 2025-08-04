### scz_gs_score.sh ###############################################################################
# specify common file paths:
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"
magma_step1_out="/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1"
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
cd /u/scratch/l/lixinzhe/revision_scratch/scz_gwas/
gunzip /u/scratch/l/lixinzhe/revision_scratch/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz
tail -n +74 /u/scratch/l/lixinzhe/revision_scratch/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv > SCZ-ALL-sumstat.tsv

# formatting:
awk 'BEGIN{ FS = OFS = "\t" } {print $2"\t"$11"\t"$14+$15}' "SCZ-ALL-sumstat.tsv" > "SCZ-ALL-sumstat_N.tsv"
new_header='SNP\tP\tN'
sed -i "1 s/.*/$new_header/" SCZ-ALL-sumstat_N.tsv

# run MAGMA to get scores:
mkdir -p /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz
${magma_dir}/magma \
    --bfile ${magma_dir}/aux/g1000_eur \
    --pval /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/SCZ-EUR-sumstat_N.tsv use='SNP,P' ncol='N' \
    --gene-annot /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/gs/out/step1.genes.annot \
    --out /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz

Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-magma-stats.R" \
    /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz.genes.out \
    ${magma_dir}/aux/NCBI37.3.gene.loc \
    "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-input.txt"

# change the traits header:
new_header="gene\tTrubetskoy_pardinas_2022"
sed -i "1 s/.*/$new_header/" "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-input.txt"

# deploy munge gs for the SCZ magma output
scdrs munge-gs \
    --out-file /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-output.gs \
    --zscore-file "/u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-input.txt" \
    --weight zscore \
    --n-max 1000 \
    --n_min 50

less /u/home/l/lixinzhe/project-geschwind/data/SCZ-GWAS/ALL/out/step2/scz-munge-output.gs
less "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" | grep PASS_Schizophrenia_Pardinas2018