### rowSum mch ####################################################################################
# load in libraries: 
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime

# load in the dataset:
data_path = '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.h5ad'
adata = sc.read_h5ad(data_path)
meta = pd.read_csv('/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_data.csv')
system_date = datetime.date.today()

### PROCESS #######################################################################################
# compute the row averge:
row_sum = np.sum(adata.X, axis = 1)
column_sum = np.sum(adata.X, axis = 0)

# subset the metas to the same order of cells:
meta.index = meta.cell
meta = meta.loc[adata.obs_names,:]

# check if the cell names are correct:
assert (adata.obs_names == meta.cell).all()

# put the columns into meta:
meta['rowSum'] = row_sum
meta['log_scale_rowsum'] = np.log1p(row_sum) - np.mean(np.log1p(row_sum))
meta['rowSum_centered'] = row_sum - np.mean(row_sum)

# also outpuut cov:
cov = pd.DataFrame({
    'cell' : meta.cell,
    'const': 1,
    'log_scale_rowsum': meta.log_scale_rowsum
})
cov.index = cov.cell

# output the cov:
cov.to_csv('/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.cov', index = False, sep = '\t')

# output the meta:
meta.to_csv('/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.csv', index = False)