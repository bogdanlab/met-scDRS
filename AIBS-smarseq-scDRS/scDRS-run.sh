out_dir="/u/project/geschwind/lixinzhe/scDRS-output/scDRS-output/AIBS-psych-trait-scDRS/with_cov/"
gs_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/"
scdrs compute-score \
    --h5ad-file "/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/expr.h5ad" \
    --h5ad-species human \
    --gs-file "${gs_dir}subset-magma_10kb_top1000_zscore.74_traits.rv1.gs" \
    --gs-species human \
    --out-folder ${out_dir} \
    --cov-file "/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/expr_covariate.cov" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
