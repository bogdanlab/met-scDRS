###########################################################################################
######                              extract methylation                              ######
###########################################################################################
# load in files:
tien_10k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien10k/Tien10k_gene_body.mcds"
tien_20k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien20k/Tien20k_gene_body.mcds"
tien_30k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien30k/Tien30k_gene_body.mcds"
tien_40k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien40k/Tien40k_gene_body.mcds"
tien_50k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien50k/Tien50k_gene_body.mcds"

# load in libraries:
import numpy as np
import pandas as pd
import xarray as xr
import anndata as ad
from scipy import sparse
from ALLCools.mcds import MCDS
from tqdm.auto import tqdm
import re
import os
import gc

# load in the mcds:
files_to_loop = [tien_10k, tien_20k, tien_30k, tien_40k, tien_50k]
output_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/'
os.makedirs(output_dir, exist_ok=True)

def get_fraction(mcsds_path, var_dim, mc_type, return_cov:bool = False):
    """
    Convert All Cools MCDS region file into a h5ad
    
    Parameters
    ----------
    mcsds_path: str
        a path to the mcsds file generated using allcools software
    var_dim: str
        feature variable, should be promoter, exon or intron
    mc_type: str
        methylation type to extract, should be either CHN or CGN
    return_cov: bool
        should the cell by feature coverage matrix be returned
    
    
    Output:
    -------
    out_h5ad:
        a h5ad file that has fraction as .X, gene names are variable, and cell names as row
        contains two meta data: coverage for each cell and gene, and count for each cell and gene
    """
    mcds = MCDS.open(
        mcsds_path,
        var_dim=var_dim,
        chunks="auto",
    )
    
    # get the data:
    da = mcds[f"{var_dim}_da"]
    mc = da.sel(mc_type=mc_type, count_type="mc")
    cov = da.sel(mc_type=mc_type, count_type="cov")
    
    # check data integrity:
    assert mc.dims == cov.dims, f"mc dims {mc.dims} != cov dims {cov.dims}"
    assert mc.sizes == cov.sizes, f"mc sizes {mc.sizes} != cov sizes {cov.sizes}"

    for dim in mc.dims:
        np.testing.assert_array_equal(
            mc[dim].values,
            cov[dim].values,
            err_msg=f"mc and cov coordinates differ on dimension {dim}",
        )
    
    if (var_dim == 'exon') or (var_dim == 'intron'):
        feature_names = pd.Index(mcds[var_dim].values)
        gene_names = feature_names.str.replace(f"_{var_dim}.*$", "", regex=True)
        unique_genes = pd.Index(gene_names).unique()
        gene_codes = pd.Categorical(gene_names, categories=unique_genes).codes
        
        n_cells = mcds.sizes["cell"]
        n_genes = len(unique_genes)
        
        mc_gene = np.zeros((n_cells, n_genes), dtype=np.float32)
        cov_gene = np.zeros((n_cells, n_genes), dtype=np.float32)
    
        batch_size = 100
        
        for start in tqdm(range(0, n_cells, batch_size), desc=f"Aggregating {var_dim} {mc_type}"):
            end = min(start + batch_size, n_cells)
            mc_block = mc.isel(cell=slice(start, end)).compute().values.astype(np.float32)
            cov_block = cov.isel(cell=slice(start, end)).compute().values.astype(np.float32)
            
            # Add interval counts into gene columns
            row_idx = np.arange(end - start)[:, None]
            np.add.at(
                mc_gene[start:end, :],
                (row_idx, gene_codes[None, :]),
                mc_block,
            )
            np.add.at(
                cov_gene[start:end, :],
                (row_idx, gene_codes[None, :]),
                cov_block,
            )
        
        frac_np = np.divide(
            mc_gene,
            cov_gene,
            out=np.zeros_like(mc_gene, dtype=np.float32),
            where=cov_gene > 0,
        )
        
        # get the meta data:
        cell_total_cov = cov_gene.sum(axis=1)
        cell_total_mc = mc_gene.sum(axis=1)
        gene_total_cov = cov_gene.sum(axis=0)
        gene_total_mc = mc_gene.sum(axis=0)
        
        # internal marker:
        var_names = unique_genes
        collapsed_from_intervals = True
        
        # set the coverage:
        cov_np = cov_gene
    
    else:
        # define the fraction:
        frac = mc / cov
        frac = frac.where(cov > 0)

        # Metadata summaries
        cell_total_cov = cov.sum(dim=var_dim).compute().values
        cell_total_mc = mc.sum(dim=var_dim).compute().values
        gene_total_cov = cov.sum(dim="cell").compute().values
        gene_total_mc = mc.sum(dim="cell").compute().values
        
        # Compute fraction matrix only after summaries
        frac_np = frac.compute().values.astype("float32")
        
        # get the variable names for naming the dataframe later:
        var_names = mcds[var_dim].values
        collapsed_from_intervals = False
        cov_np = cov.compute().values
    
    # observation
    obs = pd.DataFrame(index=mcds["cell"].values)
    obs["total_cov"] = cell_total_cov.astype("float64")
    obs["total_mc"] = cell_total_mc.astype("float64")
    obs["global_frac"] = np.divide(
        obs["total_mc"].values,
        obs["total_cov"].values,
        out=np.full(obs.shape[0], np.nan, dtype="float64"),
        where=obs["total_cov"].values > 0,
    )
    
    # get the varaible meta data:
    var = pd.DataFrame(index=var_names)
    if not collapsed_from_intervals:
        for field in ["chrom", "start", "end"]:
            coord_name = f"{var_dim}_{field}"
            if coord_name in mcds.coords:
                var[field] = mcds[coord_name].values
    else:
        var["collapsed_from_intervals"] = True
    
    # get the meta
    var["total_cov"] = gene_total_cov.astype("float64")
    var["total_mc"] = gene_total_mc.astype("float64")
    var["global_frac"] = np.divide(
        var["total_mc"].values,
        var["total_cov"].values,
        out=np.full(var.shape[0], np.nan, dtype="float64"),
        where=var["total_cov"].values > 0,
    )
    
    var["mc_type"] = mc_type
    var["region_type"] = var_dim
    
    # make sure all the shape is correct:
    expected_shape = (obs.shape[0], var.shape[0])
    assert frac_np.shape == expected_shape, (
        f"Fraction shape {frac_np.shape} != expected {expected_shape}"
    )
    
    adata = ad.AnnData(
        X=frac_np,
        obs=obs,
        var=var,
    )
    adata.X = np.nan_to_num(adata.X, nan=0.0)

    if return_cov:
        assert cov_np.shape == expected_shape, (
            f"Coverage shape {cov_np.shape} != expected {expected_shape}"
        )
        cov_adata = ad.AnnData(
            X=cov_np,
            obs=obs.copy(),
            var=var.copy(),
        )
        return adata, cov_adata
    else:
        return adata

for file in files_to_loop:
    for mc_type in ['CHN', 'CGN']:
        feature_space = 'gene_body'
        h5ad, cov_adata = get_fraction(mcsds_path = file, var_dim = feature_space, mc_type = mc_type, return_cov = True)
        output_base = re.sub('.mcds', '', re.sub('.*/', '', file))
        h5ad.write_h5ad(f"{output_dir}/{output_base}_{mc_type}_{feature_space}_extracted.h5ad")
        cov_adata.write_h5ad(f"{output_dir}/{output_base}_{mc_type}_{feature_space}_cov_only.h5ad")
        print(f'extracted data written to {output_dir}/{output_base}_{mc_type}_{feature_space}_extracted.h5ad')
        del h5ad
        del cov_adata
        gc.collect()
