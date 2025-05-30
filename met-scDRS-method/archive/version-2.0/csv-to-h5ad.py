### csv-to-h5ad.py ################################################################################
# purpose: to create a h5ad object directly from the csv file using anndata

### PREAMBLE ######################################################################################
# load in packages:
import argparse
import pandas as pd
import anndata
import scanpy as sc

# first ask for the input and output:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert CSV to h5ad')
    parser.add_argument('fraction_path', type=str, help='Path to the input CSV file')
    parser.add_argument('meta_path', type=str, help='Path to the meta data')
    parser.add_argument('h5ad_path', type=str, help='Path to the output h5ad file')
    args = parser.parse_args()

# for testing:
# fraction_path = '/u/scratch/l/lixinzhe/tmp-file/tmp-output/processed-mch.csv'
# h5ad_path = '/u/scratch/l/lixinzhe/tmp-file/tmp-output/processed-mch.h5ad'
# meta_path = '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_data.csv'

# load in the user input:
fraction_path = args.fraction_path.strip()
meta_path = args.meta_path.strip()
h5ad_path = args.h5ad_path.strip()

# data loading:
fraction = sc.read_csv(fraction_path, first_column_names = True)
meta = pd.read_csv(meta_path, index_col = 0)

### PROCESS #######################################################################################
# reorder the meta data to match the same order as meta:
fraction_rownames = fraction.obs_names
meta = meta.reindex(fraction_rownames)
meta_rownames = meta.index

# check if the cell identity are the same: 
identity_check = fraction_rownames.equals(meta_rownames)
if identity_check:
    print("Row names are identical.")
else:
    print("Row names are not identical.")

# put the meta into the adata:
for meta_field in meta.columns:
    fraction.obs[meta_field] = meta[meta_field]

# finally write the results into a h5ad file:
fraction.write(h5ad_path)
