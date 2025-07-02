import pickle
import pandas as pd
import polars as pl

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

# write out the csv file:
print('writing csv file:')
h5ad_pl.write_csv("/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.csv")
