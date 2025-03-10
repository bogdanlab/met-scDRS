### meta-process.py ################################################################################
# purpose: output meta file that contains the embedding and cell meta

### PREAMBLE ######################################################################################
import scanpy as sc
import gc
import pandas as pd

# load in merged adata:
merged_adata = sc.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/h5ad/mch/merged.h5ad')

### PROCESSING ####################################################################################
# extract out the meta data:
metadata = merged_adata.obs

# extract out the latent embedding:
umap_embedding = merged_adata.obsm['X_umap']
tsne_embedding = merged_adata.obsm['X_tsne']
tsne_df = pd.DataFrame(tsne_embedding, columns=['tSNE_1', 'tSNE_2'])
umap_df = pd.DataFrame(umap_embedding, columns=['UMAP_1', 'UMAP_2'])

# rename the rownames before merging:
tsne_df.index = metadata.index
umap_df.index = metadata.index

# Merge the embeddings with the metadata
meta = pd.concat([metadata, tsne_df, umap_df], axis = 1)

# output the merged csv:
meta.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')
