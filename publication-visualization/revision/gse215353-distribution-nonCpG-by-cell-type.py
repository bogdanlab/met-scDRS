# purpose: draw supplementary figures, get the distribution between non-CpG post process with respect to cell type

# load in data:
# load in packages:
import pickle
import pandas as pd
import polars as pl

###########################################################################################
######                                   LOADING DATA                                ######
###########################################################################################
# load in the meta data:
meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')

# load in the actual data:
intermediate_pkl = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.pkl'
with open(intermediate_pkl, 'rb') as f:
    processed_h5ad = pickle.load(f)

# use polars to output the file:
h5ad_pl = pl.DataFrame(
    processed_h5ad.X,
    schema = processed_h5ad.var_names.to_list()
    ).with_columns([
    pl.Series('cell', processed_h5ad.obs_names.to_list())
])

# move cell column to the first column
h5ad_pl = h5ad_pl.select(['cell'] + h5ad_pl.columns[:-1])
