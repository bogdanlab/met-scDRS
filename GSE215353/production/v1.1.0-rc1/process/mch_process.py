import scanpy as sc
import gc
import pandas as pd

### PROCESSING ####################################################################################
# load in merged adata:
print('reading in h5ad file')
hgnc_path = '/u/project/geschwind/lixinzhe/data/2023-08-09-mart-grch38-gene-names-dictionary.txt'
merged_adata = sc.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mch/merged.h5ad')

# first load in the mapping file:
print('mapping gene names')
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
print('checking no missingness')
gene_names = merged_adata.var_names.tolist()
gene_names_series = pd.Series(gene_names)
mapped_gene_names = gene_names_series.map(gene_map)
empty_or_na_indices = mapped_gene_names.index[mapped_gene_names.isna() | (mapped_gene_names == "")]

# rename the subset_adata to the gene names:
print('making unique')
merged_adata.var_names = pd.Index(mapped_gene_names)
merged_adata.var_names_make_unique()

# output:
print('outputting')
merged_adata.write('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-mch-v1_0_0_rc1.h5ad')
