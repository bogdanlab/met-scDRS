# load in the meta data for each of the methods:
# contain QC and output of meta data:
import anndata as ad
from tqdm.auto import tqdm
from datetime import date
today = date.today().strftime('%m%d%y')
import pandas as pd
import importlib
import met_scdrs
import numpy as np
import os
import gc

pd.set_option('display.max_columns', None)
merged_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged/'
output_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/'

# for each of the feature, we can do this preprocessing:
os.makedirs(f"{output_dir}", exist_ok = True)

features = ['gene_body']
mc_types = ['CHN', 'CGN']

pbar = tqdm(
    total=len(features) * len(mc_types),
    desc="QC h5ad files",
)

num_cells_filtered_df = pd.DataFrame(
    index=features,
    columns=mc_types,
    dtype="Int64",
)

num_genes_filtered_df = pd.DataFrame(
    index=features,
    columns=mc_types,
    dtype="Int64",
)

for feature in features:
    for mc_type in mc_types:
        h5ad_path = f"{merged_dir}/merged_07172026_{mc_type}_{feature}_fraction_raw.h5ad"
        
        # load in the h5ad path:
        tqdm.write(f"Loading {h5ad_path}")
        current_adata = ad.read_h5ad(h5ad_path)
        
        # preprocess the h5ad file:
        if mc_type == 'CGN':
            processed_h5ad, num_cell_filtered, num_gene_filtered = met_scdrs.qc_h5ad(current_adata, min_gene_coverage = 50_000)
        else:
            processed_h5ad, num_cell_filtered, num_gene_filtered = met_scdrs.qc_h5ad(current_adata, min_gene_coverage = 100_000)
        
        # document:
        num_cells_filtered_df.loc[feature, mc_type] = num_cell_filtered
        num_genes_filtered_df.loc[feature, mc_type] = num_gene_filtered

        # get the covariate file:
        row_sum = np.asarray(
            processed_h5ad.X.sum(axis=1)
        ).ravel()
        log_scale_rowsum = np.log1p(row_sum) - np.mean(np.log1p(row_sum))
        
        cov = pd.DataFrame({
            "cell": processed_h5ad.obs.index.to_numpy(),
            "const": 1,
            "log_scale_rowsum": log_scale_rowsum,
        })
        cov.to_csv(f"{output_dir}/{mc_type}_{feature}_centered_log_rowsum.cov", index = False, sep = '\t')

        # output to a preprocessed folder:
        processed_h5ad.write_h5ad(f"{output_dir}/merged_{today}_{mc_type}_{feature}_QC.h5ad")
        
        # remove the original copy:
        del current_adata
        del processed_h5ad
        gc.collect()
        
        # add to the progress bar:
        pbar.update(1)

pbar.close()

# output the summary stat
num_cells_filtered_df.to_csv(f'/u/home/l/lixinzhe/project-geschwind/plot/{today}_QC_summary_num_cells.csv')
num_genes_filtered_df.to_csv(f'/u/home/l/lixinzhe/project-geschwind/plot/{today}_QC_summary_num_genes.csv')
