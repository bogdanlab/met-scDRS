### rowSum mcg ####################################################################################
# load in libraries: 
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime

# load in the dataset:
data_path = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-mcg-v1_0_0_rc1.h5ad'
adata = sc.read_h5ad(data_path)
meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')
system_date = datetime.date.today()
subset_meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv')

### PROCESS #######################################################################################
# compute the row averge:
row_sum = np.sum(adata.X, axis = 1)
column_sum = np.sum(adata.X, axis = 0)

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

# also create one for 75k subset:
cov.index = cov.cell
subset_meta.index = subset_meta.cell
subset_cov = cov.loc[subset_meta.cell, :]

# assertion:
assert (subset_cov.index == subset_meta.index).all()

# output the cov:
cov.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/full-mcg-centered-log-library.cov', index = False, sep = '\t')
subset_cov.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/subset75k-mcg-centered-log-library.cov', index = False, sep = '\t')

# output the meta:
meta.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data_with_rowSum_mcg.csv', index = False)
