### data-extraction.py ############################################################################
# purpose: merge the seperate h5ad object, and do our met-scDRS processing

### PREAMBLE ######################################################################################
# load in packages:
import scanpy as sc
import gc
import pandas as pd

# specify the data path:
excitatory = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mcg/74eacfbe-bec4-4348-b9b4-be17630baab8.h5ad'
inhibitory = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mcg/93ea50ca-85ac-48fb-bc1f-634c247d6d64.h5ad'
non_neuronal = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mcg/6bc0a805-a1b2-4c0d-8891-cb5c8804a927.h5ad'
hgnc_path = '/u/project/geschwind/lixinzhe/data/2023-08-09-mart-grch38-gene-names-dictionary.txt'

### DATA MERGE ####################################################################################
# Load the first h5ad file
merged_adata = sc.read_h5ad(non_neuronal)

# Concatenate the next files one by one
for file in [excitatory, inhibitory]:
    print("loading in data:", file)
    adata = sc.read_h5ad(file)
    merged_adata = merged_adata.concatenate(adata, join='inner', index_unique=None)
    # Free up memory
    del adata
    gc.collect()

merged_adata.write('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mcg/merged.h5ad')

### PROCESSING ####################################################################################
# load in merged adata:
#merged_adata = sc.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mch/merged.h5ad')

# first load in the mapping file:
print('converting gene name: ')
hgnc = pd.read_csv(hgnc_path, sep = '\t')
hgnc = hgnc.drop_duplicates(subset='Gene stable ID', keep='first')

# next filter the NA and '' out:
# Remove rows where "Gene stable ID" or "Gene name" are empty strings or NA
filtered_mapping_df = hgnc[
    (hgnc['Gene stable ID'] != '') & 
    (hgnc['Gene stable ID'].notna()) & 
    (hgnc['Gene name'] != '') & 
    (hgnc['Gene name'].notna())]

# Create mapping dictionary:
gene_map = pd.Series(filtered_mapping_df['Gene name'].values, index=filtered_mapping_df['Gene stable ID'].values).to_dict()

# start mapping:
gene_names = merged_adata.var_names.tolist()
gene_names_series = pd.Series(gene_names)
mapped_gene_names = gene_names_series.map(gene_map)

# get the set of confusing genes:
empty_or_na_indices = mapped_gene_names.index[mapped_gene_names.isna() | (mapped_gene_names == "")]

# Create a mask to keep genes that are not in empty_or_na_indices
keep_genes_mask = ~merged_adata.var_names.isin(merged_adata.var_names[empty_or_na_indices])

# Apply the mask to filter out unwanted genes
merged_adata = merged_adata[:, keep_genes_mask]
gc.collect()

# now convert again and perform check to make sure there is no missingness:
gene_names = merged_adata.var_names.tolist()
gene_names_series = pd.Series(gene_names)
mapped_gene_names = gene_names_series.map(gene_map)
empty_or_na_indices = mapped_gene_names.index[mapped_gene_names.isna() | (mapped_gene_names == "")]

# rename the subset_adata to the gene names:
merged_adata.var_names = pd.Index(mapped_gene_names)
gc.collect()

# Compute the gene variances:
print('filtering low variance genes: ')
gene_variances = pd.Series(merged_adata.X.var(axis=0), index=merged_adata.var_names)
percentile_5th = gene_variances.quantile(0.05)
variance_mask = gene_variances >= percentile_5th

# filter based on the variance mask:
merged_adata = merged_adata[:, variance_mask]

# finally flip the merged_adata:
merged_adata.X = 1 - merged_adata.X

### output ########################################################################################
# write out the h5ad:
print('writing out processed h5ad object:')
merged_adata.write('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-mcg.h5ad')
