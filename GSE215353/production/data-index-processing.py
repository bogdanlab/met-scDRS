### data-index-processing.py ######################################################################
# purpose: remove duplicated index from the h5ad object:

### PREAMBLE ######################################################################################
# load in packages:
import scanpy as sc
import gc
import pandas as pd

### PROCESSING ####################################################################################
## process mcg:
# load in merged adata:
merged_adata = sc.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-mcg.h5ad')

# get the cell names:
cell_names = merged_adata.obs_names
gene_names = merged_adata.var_names

# find duplicated entries:
print(cell_names.duplicated().sum()) # return 0
print(gene_names.duplicated().sum()) # return 1193

# remove the duplicated genes:
duplicated_genes = merged_adata.var_names.duplicated()
merged_adata = merged_adata[:, ~duplicated_genes]

# output the nonduplicated reuslt:
merged_adata.write('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mcg.h5ad')

## process mch:
del merged_adata
gc.collect()

## process mch:
merged_adata = sc.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-mch.h5ad')

# get the cell names:
cell_names = merged_adata.obs_names
gene_names = merged_adata.var_names

# find duplicated entries:
print(cell_names.duplicated().sum()) # return 0
print(gene_names.duplicated().sum()) # return 1193

# remove the duplicated genes:
duplicated_genes = merged_adata.var_names.duplicated()
merged_adata = merged_adata[:, ~duplicated_genes]

# output the nonduplicated reuslt:
merged_adata.write('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mch.h5ad')
