### MC-group-test.sh ##############################################################################
# purpose: perform group testing for heterogeneity:

# first for mch:
current_date=$(date +"%Y-%m-%d")

## for cell type:
# perform downstream heterogeneity analysis for mch:
scdrs perform-downstream \
    --h5ad-file "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mch.h5ad" \
    --score-file "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/@.full_score.gz" \
    --out-folder "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/MC-downstream/GSE215353-subset/cell_type/" \
    --group-analysis X_MajorType \
    --flag-filter-data False \
    --flag-raw-count False

## for region:
# perform downstream heterogeneity analysis for mch:
scdrs perform-downstream \
    --h5ad-file "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mch.h5ad" \
    --score-file "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/@.full_score.gz" \
    --out-folder "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/MC-downstream/GSE215353-subset/region/" \
    --group-analysis X_Region \
    --flag-filter-data False \
    --flag-raw-count False
