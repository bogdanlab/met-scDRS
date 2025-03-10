# PURPOSE: try out scDRS with exted gene set on ASD expression data:
module load gcc
module load intel

### MDD GWAS MAGMA ################################################################################
# specify variables:
magma_dir="/u/project/geschwind/lixinzhe/magma/"
out_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/"
MDD_out_dir="/u/project/geschwind/lixinzhe/scDRS-output/MDD/"
trait="MDD-GWAS"
GWAS="/u/project/geschwind/lixinzhe/data/MDD-GWAS/"

# deploy scDRS score computation:
scdrs compute-score \
    --h5ad-file "${out_dir}MDD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file ${out_dir}${trait}-grch37-munge-gs-output.gs \
    --gs-species human \
    --out-folder ${MDD_out_dir} \
    --cov-file "${out_dir}MDD-covariates.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
# Trait=MDD-GWAS, n_gene=870: 0/50899 FDR<0.1 cells, 11/50899 FDR<0.2 cells (sys_time=896.9s)

# deploy plotting script:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/MDD-scDRS/scDRS-score-MDD-visualization-pipeline.R \
    "${MDD_out_dir}MDD-GWAS.score.gz" \
    "MDD-rna-GWAS" \
    "/u/project/pasaniuc/lixinzhe/plot/"

### MDD GWAS KC MAGMA #############################################################################
out_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/"
MDD_out_dir="/u/project/geschwind/lixinzhe/scDRS-output/MDD/"
gs_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/"
scdrs compute-score \
    --h5ad-file "${out_dir}MDD-RNA.h5ad" \
    --h5ad-species human \
    --gs-file "${gs_dir}subset-magma_10kb_top1000_zscore.74_traits.rv1.gs" \
    --gs-species human \
    --out-folder ${MDD_out_dir} \
    --cov-file "${out_dir}MDD-covariates.txt" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# Trait=PASS_Alzheimers_Jansen2019, n_gene=864: 5/50899 FDR<0.1 cells, 37/50899 FDR<0.2 cells (sys_time=684.2s)
# Trait=PASS_BIP_Mullins2021, n_gene=906: 7/50899 FDR<0.1 cells, 42/50899 FDR<0.2 cells (sys_time=1312.5s)
# Trait=PASS_MDD_Howard2019, n_gene=901: 0/50899 FDR<0.1 cells, 2/50899 FDR<0.2 cells (sys_time=1952.6s)
# Trait=PASS_Schizophrenia_Pardinas2018, n_gene=918: 0/50899 FDR<0.1 cells, 0/50899 FDR<0.2 cells (sys_time=2619.3s)
# Trait=UKB_460K.body_HEIGHTz, n_gene=886: 0/50899 FDR<0.1 cells, 0/50899 FDR<0.2 cells (sys_time=3188.7s)

# deploy plotting script on the Kangcheng's published gene set:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/MDD-scDRS/scDRS-score-MDD-visualization-pipeline.R \
    "${MDD_out_dir}PASS_MDD_Howard2019.full_score.gz" \
    "PASS-MDD-Howard2019-rna-GWAS" \
    "/u/project/pasaniuc/lixinzhe/plot/"
