### rowSum mcg ####################################################################################
# load in libraries: 
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime

# load in the dataset:
data_path = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mcg.h5ad'
adata = sc.read_h5ad(data_path)
meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')
system_date = datetime.date.today()

### PROCESS #######################################################################################
# compute the row averge:
row_sum = np.sum(adata.X, axis = 1)
column_sum = np.sum(adata.X, axis = 0)

# put the columns into meta:
meta['rowSum'] = row_sum

# output the meta:
meta.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data_with_rowSum_mcg.csv', index = False)