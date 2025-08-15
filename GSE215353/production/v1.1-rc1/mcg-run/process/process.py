# PURPOSE:
# output a preprocessed csv both 75K subset and full dataset for production

# load in packages:
import pickle
import pandas as pd
import polars as pl

###########################################################################################
######                                   LOADING DATA                                ######
###########################################################################################
# load in the meta data:
subset_meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv')

# load in the actual data:
intermediate_pkl = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mcg-v1_1_1_rc1.pkl'
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
h5ad_pl.write_csv("/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mcg-v1_1_1_rc1.csv")

# create an indexer from meta for joining with pl:
indexer = pl.DataFrame({
    "cell": subset_meta['cell'],
    "order": range(len(subset_meta['cell']))
})

h5ad_pl = (
    indexer
    .join(h5ad_pl, on = 'cell', how = 'left')
    .sort('order')
    .drop('order')
)

# data checking:
schema = h5ad_pl.schema

# Count float32
num_float32 = sum(dtype == pl.Float32 for dtype in schema.values())

# Count strings
num_str = sum(dtype == pl.Utf8 for dtype in schema.values())

# assertion:
assert num_str == 1, 'more than one column of str detected!'
assert ((len(schema) - 1) == num_float32), 'not all other columns are float32!'

# write the file:
h5ad_pl.write_csv("/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-75K-subset-mcg-v1_1_1_rc1.csv")
