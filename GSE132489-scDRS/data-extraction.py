"""

Usage:
    data-extraction.py [--input_dir=<in_dir> --output_dir=<out_dir>]

Options:
    --input_dir=<in_dir>      input directory for a list of hdf5 files
    --output_dir=<out_dir>    output directory for concatenated files

"""
from docopt import docopt
import xarray as xr
import pandas as pd
import glob
import os

# define paths:
arguments = docopt(__doc__)
print(arguments)
data_dir = arguments['--input_dir']
out_dir = arguments['--output_dir']

### LOAD DATA #####################################################################################
# load hdf5:
hdf5_file = glob.glob(data_dir + '*.hdf5')
mcds = xr.open_mfdataset(hdf5_file, concat_dim = 'cell', combine = 'nested')

# get mch 100kb data:
chrom100k_data = mcds['chrom100k_da']

# for mch:
mch_count = chrom100k_data.sel(mc_type = 'CHN', count_type = 'mc').squeeze()
mch_coverage = chrom100k_data.sel(mc_type = 'CHN', count_type = 'cov').squeeze()

# for mcg:
mcg_count = chrom100k_data.sel(mc_type = 'CGN', count_type = 'mc').squeeze()
mcg_coverage = chrom100k_data.sel(mc_type = 'CGN', count_type = 'cov').squeeze()

# to pandas:
mcg_count_df = mcg_count.to_pandas()
mch_count_df = mch_count.to_pandas()
mcg_coverage_df = mcg_coverage.to_pandas()
mch_coverage_df = mch_coverage.to_pandas()

### GENE LEVEL ####################################################################################
gene_data = mcds['gene_da']
# get mcg:
mcg_gene_level_count = gene_data.sel(mc_type = 'CGN', count_type = 'mc').squeeze()
mcg_gene_level_coverage = gene_data.sel(mc_type = 'CGN', count_type = 'cov').squeeze()

# get mch:
mch_gene_level_count = gene_data.sel(mc_type = 'CHN', count_type = 'mc').squeeze()
mch_gene_level_coverage = gene_data.sel(mc_type = 'CHN', count_type = 'cov').squeeze()

# convert to pandas:
mcg_gene_level_count_df = mcg_gene_level_count.to_pandas()
mcg_gene_level_coverage_df = mcg_gene_level_coverage.to_pandas()
mch_gene_level_count_df = mch_gene_level_count.to_pandas()
mch_gene_level_coverage_df = mch_gene_level_coverage.to_pandas()

# OUTPUT:
# 100kb:
mch_count_df.to_csv(out_dir + 'full_100_kb_mch_count' + '.csv')
mch_coverage_df.to_csv(out_dir + 'full_100_kb_mch_coverage' + '.csv')
mcg_count_df.to_csv(out_dir + 'full_100_kb_mcg_count' + '.csv')
mcg_coverage_df.to_csv(out_dir + 'full_100_kb_mcg_coverage' + '.csv')

# gene level:
mch_gene_level_count_df.to_csv(out_dir + 'full_gene_mch_count' + '.csv')
mch_gene_level_coverage_df.to_csv(out_dir + 'full_gene_mch_coverage' + '.csv')
mcg_gene_level_count_df.to_csv(out_dir + 'full_gene_mcg_count' + '.csv')
mcg_gene_level_coverage_df.to_csv(out_dir + 'full_gene_mcg_coverage' + '.csv')
