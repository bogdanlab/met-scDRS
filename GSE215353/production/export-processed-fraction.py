### export-processed-fraction.py ##################################################################
# purpose: export the processed fraction from the post-processed h5ad file to a csv file:

### PREAMBLE ######################################################################################
# load in libraries:
import anndata
import pandas as pd
import gc

# load in the file:
mch = anndata.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mch.h5ad')

# Convert the data into pandas dataframe:
fraction_df = pd.DataFrame(mch.X, index = mch.obs_names, columns = mch.var_names)

# write the data out to csv:
fraction_df.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mch.csv')

# clear memory:
del mch
del fraction_df
gc.collect()

# load in the file:
mcg = anndata.read_h5ad('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mcg.h5ad')

# Convert the data into pandas dataframe:
fraction_df = pd.DataFrame(mcg.X, index = mcg.obs_names, columns = mcg.var_names)

# write the data out to csv:
fraction_df.to_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mcg.csv')
