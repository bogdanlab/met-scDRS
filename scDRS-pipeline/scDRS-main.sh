# qrsh -l h_data=8G,h_rt=72:00:00,highp -pe shared 8
module load gcc
module load intel

### PART 1: COMPUTE LOCATION FILE USING MAGMA #####################################################
# Step 1: download MAGMA software, gene location file, and reference data from 
# https://ctg.cncr.nl/software/magma after this step, one should have a folder <MAGMA_DIR>
# with the following files:
# 1) <MAGMA_DIR>/magma 2) <MAGMA_DIR>/g1000_eur.(bed|bim|fam) 3) <MAGMA_DIR>/NCBI37.3.gene.loc

# specify variables:
magma_dir="/u/project/geschwind/lixinzhe/magma/"
out_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/"

# Step 2: make gene annotation file for MAGMA using the following command, this only needs to be done
# once for different GWAS summary statistics, and the file will be saved to out/step1.genes.annot
${magma_dir}software/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}auxiliary/ref/eur/g1000_eur.bim \
    --gene-loc ${magma_dir}auxiliary/grch37/NCBI37.3.gene.loc \
    --out ${out_dir}grch37_10kb_eur

# Step 3: run MAGMA using the following command, this takes a GWAS file ${trait}.pval,
# which at least has the columns: SNP, P, N, which corresponds to the SNP id
# (matched to the ${magma_dir}/g1000_eur.bim), p-value, sample size. For example,
# <trait>.pval file looks like
#
# CHR     BP      SNP             P           N
# 1       717587  rs144155419     0.453345    279949
# 1       740284  rs61770167      0.921906    282079
# 1       769223  rs60320384      0.059349    281744
#
# After this step, one should obtain a file out/step2/${trait}.gene.out, and the top genes with
# largest Z-scores can be input to scDRS.

### PART 2: COMPUTE GENESET FROM CUSTOM GWAS: #####################################################
# specify variables:
trait="iPSYCH-PGC_ASD_Nov2017"
GWAS="/u/project/geschwind/lixinzhe/data/ASD-GWAS/"

# just use N = # control + # patients for the GWAS for all snps:
awk 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "N" : "46350") }' "${trait}" > "${trait}_N"

# call magma to compute gene set z score and p values using custom GWAS:
${magma_dir}software/magma \
    --bfile ${magma_dir}auxiliary/ref/eur/g1000_eur \
    --pval ${GWAS}${trait}_N use='SNP,P' ncol='N' \
    --gene-annot ${out_dir}grch37_10kb_eur.genes.annot \
    --out ${out_dir}${trait}-grch37

# prepare the z-score file:
Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-magma-stats.R" \
    "${out_dir}${trait}-grch37.genes.out" \
    "${magma_dir}auxiliary/grch37/NCBI37.3.gene.loc" \
    "${out_dir}${trait}-grch37-munge-input.txt"

# deploy gs munge:
scdrs munge-gs \
    --out-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --zscore-file ${out_dir}${trait}-grch37-munge-input.txt \
    --weight zscore \
    --n-max 1000

### PART 3: COMPUTE scDRS SCORE: ##################################################################
# conversion from seurat to h5ad:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/seurat-h5ad.R \
    "/u/project/geschwind/lixinzhe/data/2023-02-17-combined-obj-preprocessed.rds" \
    "${out_dir}ASD-RNA.h5ad"

# get covariate file:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-covariates.R \
    "/u/project/geschwind/lixinzhe/data/2023-02-17-combined-obj-preprocessed.rds" \
    "~ 0 + technician + SequencingPool + ancestry + BrainBank + filled.pmi + filled.rin + Age + percent.mt" \
    "${out_dir}ASD-covariates.txt"

# compute scdrs score:
ASD_out_dir="/u/project/geschwind/lixinzhe/scDRS-output/ASD/"
scdrs compute-score \
    --h5ad-file "${out_dir}ASD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --gs-species human \
    --out-folder ${ASD_out_dir} \
    --cov-file "${out_dir}ASD-covariates.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# Trait=trait, n_gene=793: 0/18666 FDR<0.1 cells, 0/18666 FDR<0.2 cells (sys_time=505.6s)

### PART 4: Results Visualization #################################################################
# plot out the p values distribution
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/plot-pvalue-distribution.R \
    ${ASD_out_dir}"trait.score.gz" \
    "/u/project/geschwind/lixinzhe/data/2023-04-11-meta-data.txt" \
    "Diagnosis" \
    "/u/project/geschwind/lixinzhe/plot/2023-05-25-scDRS-score-distribution.png"

Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/ASD-bimodal-scDRS/scDRS-score-visualization.R" \
    "${ASD_out_dir}trait.score.gz" \
    "ASD-GWAS" \
    "/u/project/pasaniuc/lixinzhe/plot/"

### testing if the result will change if I have index as my first column's name:new_header="gene\t${trait}"
new_header="index\tconst\ttechnician\tSequencingPoolYZCL16\tSequencingPoolYZCL17\tSequencingPoolYZCL19\tSequencingPoolYZCL28\tancestryAMR\tancestryEUR\tBrainBankUM-BTB\tfilled.pmi\tfilled.rin\tAge\tpercent.mt"
sed "1 s/.*/$new_header/" "${out_dir}ASD-covariates.txt" > "${out_dir}-index-ASD-covariates.txt"

scdrs compute-score \
    --h5ad-file "${out_dir}ASD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --gs-species human \
    --out-folder ${ASD_out_dir}cov-index/ \
    --cov-file "${out_dir}-index-ASD-covariates.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
# Trait=trait, n_gene=793: 0/18666 FDR<0.1 cells, 0/18666 FDR<0.2 cells (sys_time=518.3s)

# the two results are identical
gunzip -c ${ASD_out_dir}cov-index/trait.score.gz > ${ASD_out_dir}cov-index/trait.score
gunzip -c ${ASD_out_dir}trait.score.gz > ${ASD_out_dir}trait.score
cmp ${ASD_out_dir}trait.score ${ASD_out_dir}cov-index/trait.score

### RUN ANALYSIS WITH FULL SET OF COVARIATES ######################################################
result_dir="${ASD_out_dir}full-cov/"
trait="iPSYCH-PGC_ASD_Nov2017"

scdrs compute-score \
    --h5ad-file "${out_dir}ASD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --gs-species human \
    --out-folder ${result_dir} \
    --cov-file "${out_dir}full-expression-ASD-covariates.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# Trait=trait, n_gene=793: 0/18666 FDR<0.1 cells, 0/18666 FDR<0.2 cells (sys_time=495.2s)

# visualize result:
Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/ASD-bimodal-scDRS/scDRS-score-visualization.R" \
    "${result_dir}trait.score.gz" \
    "ASD-GWAS-full-covariate" \
    "/u/project/pasaniuc/lixinzhe/plot/"

### ADD GENE NUMBER COVARIATES ####################################################################
result_dir="${ASD_out_dir}nfeature-cov/"
trait="iPSYCH-PGC_ASD_Nov2017"

scdrs compute-score \
    --h5ad-file "${out_dir}ASD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --gs-species human \
    --out-folder ${result_dir} \
    --cov-file "${out_dir}ASD-covariates-with-nFeature.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# Trait=trait, n_gene=793: 0/18666 FDR<0.1 cells, 0/18666 FDR<0.2 cells (sys_time=500.6s)

#visualize result:
Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/ASD-bimodal-scDRS/scDRS-score-visualization.R" \
    "${result_dir}trait.score.gz" \
    "ASD-GWAS-with-nfeature-covariate" \
    "/u/project/pasaniuc/lixinzhe/plot/"
