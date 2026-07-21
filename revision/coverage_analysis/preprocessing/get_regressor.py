# obtain regressor for regressing cell level batch effect
# load in packages:
import anndata as ad
import re
import os
import gc
import pandas as pd
from tqdm.auto import tqdm
import numpy as np
from datetime import date
import met_scdrs
import importlib
today = date.today().strftime("%m%d%Y")
from statsmodels.stats.outliers_influence import variance_inflation_factor

###########################################################################################
######                                  Define function                              ######
###########################################################################################
def generate_control_metric(
    frac_adata,
    coverage_adata,
    covariate_file_output: str,
    total_coverage: str = 'total_cov',
    total_mc: str = 'total_mc',
    ):
    """
    Purpose: generate control meta data such as missing rate, coverage density, cell coverage, cell read depth
    
    Parameters
    ----------
    frac_adata: anndata
        adata.X contains the fraction of methylation per cell per gene
        adata.obs contains cell level meta data
        adata.var contains gene level meta data
    coverage_adata: anndata
        coverage_adata.X contains the raw coverage count of methylation per cell per gene
    
    """
    
    # check the inputt, make sure there is no missing gene and cell from the coverage
    missing_cells = frac_adata.obs_names.difference(coverage_adata.obs_names)
    missing_genes = frac_adata.var_names.difference(coverage_adata.var_names)

    if len(missing_cells) > 0 or len(missing_genes) > 0:
        raise ValueError(
            f"Missing {len(missing_cells)} cells and {len(missing_genes)} genes "
            "from coverage_adata."
        )

    coverage_adata = coverage_adata[
        frac_adata.obs_names,
        frac_adata.var_names
    ].copy()
    
    # calculate missiness:
    missingness = np.asarray(np.sum(coverage_adata.X == 0, axis = 1)).ravel()
    
    # prepare the covariate output matrix:
    covariate_column = ['cell', 'const', 'coverage', 'met_count', 'missing_vals']
    out_df = pd.DataFrame(index = frac_adata.obs_names, columns = covariate_column)
    
    # fill in the covariate matrix:
    out_df.loc[frac_adata.obs_names, 'cell'] = frac_adata.obs_names
    out_df.loc[frac_adata.obs_names, 'const'] = 1
    out_df.loc[frac_adata.obs_names, 'coverage'] = np.log1p(frac_adata.obs.loc[:, total_coverage])
    out_df.loc[frac_adata.obs_names, 'met_count'] = np.log1p(frac_adata.obs.loc[:, total_mc])
    out_df.loc[frac_adata.obs_names, 'missing_vals'] = np.log1p(missingness)
    
    out_df["coverage"] = out_df["coverage"].astype(float)
    out_df["met_count"] = out_df["met_count"].astype(float)
    count_coverage_correlation = np.corrcoef(out_df.coverage.values, out_df.met_count.values)[0, 1]
    print(f'\ncorrelation coef between log count and log coverage {count_coverage_correlation:.2f}')

    # check the Variance inflation factor:
    X = out_df[["const", "coverage", "met_count", "missing_vals"]].astype(float)
    vif = pd.DataFrame({
        "variable": X.columns,
        "VIF": [
            variance_inflation_factor(X.values, i)
            for i in range(X.shape[1])
        ]
    })
    print(vif)

    # output the dataframe:
    out_df.to_csv(covariate_file_output, index = False, sep = '\t')
    return out_df

###########################################################################################
######                                    non-CpG                                    ######
###########################################################################################
# load in the coverage and fraction:
coverage = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged/merged_07172026_CHN_gene_body_coverage_raw.h5ad')
#fraction = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged/merged_07172026_CHN_gene_body_fraction_raw.h5ad')
qced_fraction = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CHN_gene_body_QC.h5ad')

# output the out_df:
cov_df = generate_control_metric(
    qced_fraction,
    coverage,
    covariate_file_output = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_full.cov'
    )

###########################################################################################
######                                    CpG                                        ######
###########################################################################################
coverage = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged/merged_07172026_CGN_gene_body_coverage_raw.h5ad')
#fraction = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged/merged_07172026_CGN_gene_body_fraction_raw.h5ad')
qced_fraction = ad.read_h5ad('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/merged_071726_CGN_gene_body_QC.h5ad')

# output the out_df:
cov_df = generate_control_metric(
    qced_fraction,
    coverage,
    covariate_file_output = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CGN_gene_body_full.cov'
    )

cov_df = cov_df.drop(columns = ['met_count'])
cov_df.to_csv('/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CGN_gene_body_full.cov', index = False, sep = '\t')

