executable="/u/home/l/lixinzhe/project-github/met-scDRS/GSE215353/revision/v1.1-rc-subset/generate_data/csv_to_h5ad.py"
subset_file="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v1.2/processed-mcg-gene-name.csv"
meta_file="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv"
h5ad_file="/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v1.2/processed-75K-subset-mcg.h5ad"

python ${executable} \
    ${subset_file} \
    ${meta_file} \
    ${h5ad_file}