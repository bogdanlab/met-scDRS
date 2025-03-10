### fraction-scdrs.sh #############################################################################
# PURPOSE: try out scDRS with exted gene set on ASD expression data:
module load gcc
module load intel

# specify variables:
magma_dir="/u/project/geschwind/lixinzhe/magma/"
input_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/"
methylation_out_dir="/u/project/geschwind/lixinzhe/scDRS-output/fraction-gse132489-methyl/"

# specify scdrs run:
gs_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/"
scdrs compute-score \
    --h5ad-file "${input_dir}subset-GSE132489-mcg.h5ad" \
    --h5ad-species mouse \
    --gs-file "${gs_dir}magma_10kb_top1000_zscore.74_traits.rv1.gs" \
    --gs-species human \
    --out-folder "${methylation_out_dir}/mcg/KC/all/" \
    --cov-file "${input_dir}subset-GSE132489-covariates-mcg.txt" \
    --flag-filter-data False \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
