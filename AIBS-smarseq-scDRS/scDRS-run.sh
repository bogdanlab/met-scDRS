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

out_dir="/u/project/geschwind/lixinzhe/scDRS-output/scDRS-output/AIBS-psych-trait-scDRS/without_cov/"
gs_dir="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/"
mkdir -p $out_dir
scdrs compute-score \
    --h5ad-file "/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/expr.h5ad" \
    --h5ad-species human \
    --gs-file "${gs_dir}subset-magma_10kb_top1000_zscore.74_traits.rv1.gs" \
    --gs-species human \
    --out-folder ${out_dir} \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
# Trait=PASS_Alzheimers_Jansen2019, n_gene=966: 540/47432 FDR<0.1 cells, 1054/47432 FDR<0.2 cells (sys_time=1297.9s)
# Computing control scores: 100%|█| 1000/1000 [17:46<00:00,  1.07s/it
# Trait=PASS_BIP_Mullins2021, n_gene=968: 4/47432 FDR<0.1 cells, 157/47432 FDR<0.2 cells (sys_time=2510.0s)
# Computing control scores: 100%|█| 1000/1000 [17:57<00:00,  1.08s/it
# Trait=PASS_MDD_Howard2019, n_gene=978: 0/47432 FDR<0.1 cells, 0/47432 FDR<0.2 cells (sys_time=3734.1s)
# Computing control scores: 100%|█| 1000/1000 [18:14<00:00,  1.09s/it
# Trait=PASS_Schizophrenia_Pardinas2018, n_gene=977: 14/47432 FDR<0.1 cells, 30/47432 FDR<0.2 cells (sys_time=4972.5s)
# Computing control scores: 100%|█| 1000/1000 [17:00<00:00,  1.02s/it
# Trait=UKB_460K.body_HEIGHTz, n_gene=983: 88/47432 FDR<0.1 cells, 601/47432 FDR<0.2 cells (sys_time=6137.4s)
