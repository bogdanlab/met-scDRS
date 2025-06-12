import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from skmisc.loess import loess
from typing import List
import warnings
import time
from met_scdrs.util import get_memory
import met_scdrs
from anndata import AnnData
import gc
from tqdm import tqdm
import os

def normalize(
    h5ad_obj,
    method : str = 'inverse',
    variance_clip : int = 5,
    transformation : str = None,
    verbose = True
):
    """
    Preprocess the methylation cell by gene data
    inverse the fraction to 1 - fraction such that higher methylation corresponds to higher expected gene expression
    
    Normalization schematics:
        logit: logit transformation of proportion logit(proportion)
        arcsine transformation: arcsine(sqrt(proportion))
        library size normalization: log(1+proportion / sum(proportion))
    
    Parameters
    ----------
    h5ad_obj : anndata.AnnData
        AnnData object loaded from h5ad_file path
    method : str
        methodology to preprocess the single cell methylation matrix, supported methods:
        "inverse" : inverse the fraction into 1 - X
    transformation : str, optional
        "logit" : logit transformation
        "arcsine" : arcsine transformation
        "log_library" : inverse the fraction then library normalize
        if None, no transformation is applied
    variance_clip : int
        only genes with greater than specified percentile will be retained
        Default is 5 (i.e., remove low-variance genes under the 5th percentile).
    verbose : bool
    """
    # print header
    print(f'\n\n\n PREPROCESSING')
    initial_time = time.time()
    header = "--preprocess_method %s \\\n" % method
    header += "--variance_clip %s \\\n" % variance_clip
    header += "--transformation %s \\\n" % transforamtion
    print(header)
    
    # obtain the memory usage:
    print(f"Initiating preprocess, memory usage: {get_memory():.2f} MB") if verbose else None
    
    ###########################################################################################
    ######                                  helper function                              ######
    ###########################################################################################
    # define method (inverse): 
    def inverse_method(h5ad_obj):
        if sparse.issparse(h5ad_obj.X):
            print('h5ad object is recognized as sparse object')
            h5ad_obj.X = 1 - h5ad_obj.X.toarray()
            return h5ad_obj
        
        elif isinstance(h5ad_obj.X, np.ndarray):
            print('h5ad object is recognized as numpy dense array')
            np.subtract(1, h5ad_obj.X, out = h5ad_obj.X)
            return h5ad_obj
        
        else:
            print('detected adata.X not numpy or scipy sparse matrix')
            print('highly recommended to create h5ad with scipy sparse object, attemping 1 - h5ad_obj.X')
            h5ad_obj.X = 1 - h5ad_obj.X
            return h5ad_obj
    
    def compute_variance(h5ad_obj):
        if sparse.issparse(h5ad_obj.X):
            variances_ = np.var(h5ad_obj.X.toarray(), axis = 0)
        else:
            variances_ = np.var(h5ad_obj.X, axis = 0)
        return variances_
    
    def arcsine_transformation(adata, chunk_size):
        # in case of float point error:
        eps = 1e-8
        
        for chunk in range(0, adata.shape[0], chunk_size):
            # assert that all entries in chunk is between 0 and 1:
            assert np.all((adata.X[chunk : chunk + chunk_size] >= -eps) & (adata.X[chunk : chunk + chunk_size] <= 1 + eps))
            
            # if bound between 0 and 1, compute arcsin transformation:
            adata.X[chunk : chunk + chunk_size] = np.arcsin(np.sqrt(adata.X[chunk : chunk + chunk_size])).astype(np.float32)
    
    def logit_transformation(adata, chunk_size):
        # because 1e-8 still drive the result to infinity after np.chunk
        eps = 1e-7 
        
        # lets make a logit transformation:
        for chunk in range(0, adata.shape[0], chunk_size):
            # assert that all entries in chunk is between 0 and 1:
            assert np.all((adata.X[chunk : chunk + chunk_size] >= -eps) & (adata.X[chunk : chunk + chunk_size] <= 1 + eps))
            
            # clip current chunk:
            adata_chunk = np.clip(adata.X[chunk : chunk + chunk_size], eps, 1 - eps) # stabilize logit transformation
            adata_chunk = sp.special.logit(adata_chunk).astype(np.float32)
            
            # if bound between 0 and 1, compute logit:
            adata.X[chunk : chunk + chunk_size] = adata_chunk
    
    def library_size_normalization(adata, chunk_size):
        # in case of float point error:
        eps = 1e-8
        scale_factor = 10000
        
        # do this in chunks:
        for chunk in range(0, adata.shape[0], chunk_size):
            # assert that all entries in chunk between 0 and 1:
            assert np.all((adata.X[chunk : chunk + chunk_size] >= -eps) & (adata.X[chunk : chunk + chunk_size] <= 1 + eps))
            
            # if bound between 0 and 1, library size normalization:
            adata_chunk = adata.X[chunk : chunk + chunk_size]
            library_size = np.sum(adata_chunk, axis = 1, keepdims = True)
            adata_chunk = adata_chunk / library_size * scale_factor
            adata_chunk = np.log1p(adata_chunk).astype(np.float32)
            
            # fill in the data:
            adata.X[chunk : chunk + chunk_size] = adata_chunk
    
    ###########################################################################################
    ######                                    preprocess .X                              ######
    ###########################################################################################
    # preprocess the .X data:
    methods = {'inverse' : inverse_method}
    preprocessed_data = methods[method](h5ad_obj)
    
    gc.collect()
    
    ###########################################################################################
    ######                                    variance clip                              ######
    ###########################################################################################
    # variance masking and assertion:
    print('Filtering gene based on variance') if verbose else None
    gene_variances = compute_variance(preprocessed_data)
    var_threshold = np.percentile(gene_variances, variance_clip)
    print(f"Variance clip at {variance_clip} percentile: {var_threshold:.5f}") if verbose else None
    assert var_threshold > 0, "all gene variance is not > 0, please increase the variance clipping percentile"
    
    # filter off genes with small variance:
    high_var_mask = gene_variances >= var_threshold
    print(f"Retaining {np.sum(high_var_mask)} / {len(gene_variances)} genes")
    preprocessed_data._inplace_subset_var(high_var_mask)
    
    # assertion variance filter and gene with high variance is retained:
    preprocessed_var = compute_variance(preprocessed_data)
    assert (preprocessed_var >= var_threshold).all, "variance is not clipped properly, exitting"
    gc.collect()
    
    ###########################################################################################
    ######                                   transformation                              ######
    ###########################################################################################
    # establish transforamtion methodologies:
    transformation_method = {
        'arcsine' : arcsine_transformation,
        'logit' : logit_transformation,
        'library' : library_size_normalization
        }
    
    # if transformation is provided
    if transformation:
        assert transformation in transformation_method.keys(), 'unsupported --transformation input'
        transformation_method[transformation](preprocessed_data, chunk_size = 1000)
    
    # method to usage 
    # print our usage information:
    print(f'Normalization completed, elapsed time: {(time.time() - initial_time):.3f} seconds') if verbose else None
    print(f"Finished normalization, memory usage: {get_memory():.2f} MB") if verbose else None
    print(f"Top 5 preprocessed cells and genes: {preprocessed_data.X[0:5, 0:5]}") if verbose else None
    return preprocessed_data

def category2dummy(
    df: pd.DataFrame, cols: List[str] = None, verbose: bool = False
) -> pd.DataFrame:
    """
    Convert categorical variables in a dataframe to binary dummy variables.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to convert
    cols : List[str], optional
        Columns to convert, by default None (columns are then selected automatically)
    verbose : bool, optional
        Print progress, by default False

    Returns
    -------
    pd.DataFrame
        Converted dataframe
    """
    df = df.copy()

    if cols is None:
        # infer non-numerical columns
        cols = list(set(df.columns) - set(df._get_numeric_data().columns))

    assert set(cols).issubset(set(df.columns)), "'cols' must be a subset of df.columns"

    cols_to_drop = []
    cols_to_add = []
    dummy_dfs = []
    for col in cols:
        # create df of dummy variables
        dummy_df = pd.get_dummies(df[col], drop_first=True)
        dummy_df.columns = [f"{col}_{s}" for s in dummy_df.columns]
        dummy_df.loc[df[col].isnull(), dummy_df.columns] = np.nan
        cols_to_drop.append(col)
        cols_to_add.extend(dummy_df.columns)
        dummy_dfs.append(dummy_df)

    # all columns in dummy_dfs
    df = pd.concat([df, *dummy_dfs], axis=1).drop(columns=cols_to_drop)

    if (len(cols_to_add) > 0) and verbose:
        print(
            "scdrs.pp.category2dummy: "
            f"Detected categorical columns: {','.join(cols)}. "
            f"Added dummy columns: {','.join(cols_to_add)}. "
            f"Dropped columns: {','.join(cols_to_drop)}."
        )
    return df

def preprocess(
    data,
    weight_option,
    ctrl_match_key,
    cov=None,
    n_mean_bin=10,
    n_var_bin=10,
    n_length_bin=10,
    n_chunk=None,
    copy=False,
    verbose = True,
):
    """
    Preprocess single-cell data for scDRS analysis.

        1. Correct covariates by regressing out the covariates (including
        a constant term) and adding back the original mean for each gene.

        2. Compute gene-level and cell-level statistics for the
        covariate-corrected data.

    Information is stored in `data.uns["SCDRS_PARAM"]`. It operates in
    implicit-covariate-correction mode when `data.X` is sparse and `cov`
    not `None` to improve memory efficiency; it operates in normal mode
    otherwise.

    In normal mode, `data.X` is replaced by the covariate-corrected data.

    In implicit-covariate-correction mode, the covariate correction information
    is stored in `data.uns["SCDRS_PARAM"]` but is not explicitly applied to
    `data.X`, so that `data.X` is always sparse. Subsequent computations on
    the covariate-corrected data are based on the original data `data.X` and
    the covariate correction information. Specifically,

        CORRECTED_X = data.X + COV_MAT * COV_BETA + COV_GENE_MEAN

    Parameters
    ----------
    data : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed normalized
    weight_option : str
        weighting option, same as cli input
    ctrl_match_key : str, default="mean_var"
        Gene-level statistic used for matching control and disease genes;
        should be in `data.uns["SCDRS_PARAM"]["GENE_STATS"]`.
    cov : pandas.DataFrame, default=None
        Covariates of shape (n_cell, n_cov). Should contain
        a constant term and have values for at least 75% cells.
    n_mean_bin : int, default=10
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=10
        Number of expression-variance bins for matching control genes.
    n_length_bin : int, default = 10
        Number of gene length bins for matching control genes.
    n_chunk : int, default=None
        Number of chunks to split the data into when computing mean and variance
        using _get_mean_var_implicit_cov_corr. If n_chunk is None, set to 5/sparsity.
    copy : bool, default=False
        Return a copy instead of writing to data.
    verbose: bool, default=True


    Returns
    -------
    Overview:
        `data.X` will be updated as the covariate-corrected data in normal mode
        and will stay untouched in the implicit covariate correctoin mode.
        Preprocessing information is stored in `data.uns["SCDRS_PARAM"]`.
    FLAG_SPARSE : bool
        If data.X is sparse.
    FLAG_COV : bool
        If covariate correction is performed.
    COV_MAT : pandas.DataFrame
        Covariate matrix of shape (n_cell, n_cov).
    COV_BETA: pandas.DataFrame
        Covariate effect sizes of shape (n_gene, n_cov).
    COV_GENE_MEAN: pandas.Series
        Gene-level mean expression.
    GENE_STATS : pandas.DataFrame
        Gene-level statistics of shape (n_gene, 7):

        - "mean" : mean expression in log scale.
        - "var" : expression variance in log scale.
        - "var_tech" : technical variance in log scale. if weighting option is vs
        - "ct_mean" : mean expression in original non-log scale.
        - "ct_var" : expression variance in original non-log scale.
        - "ct_var_tech" : technical variance in original non-log scale. if weighting option is vs
        - "mean_var" : n_mean_bin * n_var_bin mean-variance bins

    CELL_STATS : pandas.DataFrame
        Cell-level statistics of shape (n_cell, 2):

        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.


    Notes
    -----
    Covariate regression:
        adata.X =  cov * beta + resid_X.
    scDRS saves:
        COV_MAT = cov, COV_BETA = (-beta), COV_GENE_MEAN = adata.X.mean(axis=0)
    The scDRS covariate-corrected data:
        CORRECTED_X = resid_X + GENE_MEAN = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN.


    """
    # print the memory usage:
    print(f"Starting preprocessing, memory usage: {get_memory():.2f} MB") if verbose else None
    tracker = met_scdrs.util.MemoryTracker()
    tracker.start()
    
    adata = data.copy() if copy else data
    n_cell, n_gene = adata.shape

    # Parameters and flags
    flag_sparse = sparse.issparse(adata.X)
    flag_cov = cov is not None
    adata.uns["SCDRS_PARAM"] = {
        "FLAG_SPARSE": flag_sparse,
        "FLAG_COV": flag_cov,
    }

    # Update adata.X
    if flag_sparse:
        # Force sparse.csr_matrix for the sparse mode
        if not isinstance(adata.X, sparse.csr_matrix):
            adata.X = sparse.csr_matrix(adata.X)
    else:
        # Force np.ndarray for the dense mode
        if not isinstance(adata.X, np.ndarray):
            adata.X = np.array(adata.X)

    # Covariate correction
    if flag_cov:
        # Check if cells in data and cov are consistent
        assert (
            len(set(cov.index) & set(adata.obs_names)) > 0.75 * n_cell
        ), "cov does not match the cells in data"

        df_cov = pd.DataFrame(index=adata.obs_names)
        df_cov = df_cov.join(cov)
        df_cov = category2dummy(df=df_cov, verbose=True)
        df_cov.fillna(df_cov.mean(), inplace=True)

        # Add const term if df_cov does not already have it (or a linear combination of it)
        v_resid = _reg_out_inplace(np.ones(n_cell), df_cov.values)
        if (v_resid ** 2).mean() > 0.01:
            df_cov["SCDRS_CONST"] = 1

        # Gene mean: numpy.ndarray of shape (n_gene,)
        v_gene_mean = np.array(adata.X.mean(axis=0)).flatten()
        if flag_sparse:
            # Sparse mode: save correction information
            mat_beta = np.linalg.solve(
                np.dot(df_cov.values.T, df_cov.values) / n_cell,
                sparse.csr_matrix.dot(df_cov.values.T, adata.X) / n_cell,
            )

            adata.uns["SCDRS_PARAM"]["COV_MAT"] = df_cov
            adata.uns["SCDRS_PARAM"]["COV_BETA"] = pd.DataFrame(
                -mat_beta.T, index=adata.var_names, columns=df_cov.columns
            )
            adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"] = pd.Series(
                v_gene_mean, index=adata.var_names
            )
        else:
            # Dense mode: regress out covariate and add back mean
            peak = tracker.stop()
            print(f'peak memory before regression: {peak:.2f} GB') if verbose else None
            
            tracker.start()
            adata.X = _reg_out_inplace(adata.X, df_cov.values)
            adata.X += v_gene_mean
            peak = tracker.stop()
            
            print(f'peak memory during regression: {peak:.2f} GB') if verbose else None
    tracker.start()
    
    
    #             # Note: this version (dense+cov) should produce the exact toydata results
    #             adata.var["mean"] = adata.X.mean(axis=0).T
    #             adata.X -= adata.var["mean"].values
    #             adata.X = reg_out(adata.X, df_cov[['SCDRS_CONST', 'covariate']].values)
    #             adata.X += adata.var["mean"]

    # Precompute for each gene and mean&var for each cell
    if flag_sparse and flag_cov:
        implicit_cov_corr = True
        if n_chunk is None:
            n_chunk = 5 * adata.shape[0] * adata.shape[1] // adata.X.data.shape[0] + 1
    else:
        implicit_cov_corr = False
        if n_chunk is None:
            n_chunk = 20
    
    cell_weight = None
    
    df_gene, df_cell = compute_stats(
        adata,
        weight_option = weight_option,
        ctrl_match_key = ctrl_match_key,
        implicit_cov_corr=implicit_cov_corr,
        cell_weight=cell_weight,
        n_mean_bin=n_mean_bin,
        n_var_bin=n_var_bin,
        n_length_bin=n_length_bin,
        n_chunk=n_chunk,
        verbose = verbose
    )

    adata.uns["SCDRS_PARAM"]["GENE_STATS"] = df_gene
    adata.uns["SCDRS_PARAM"]["CELL_STATS"] = df_cell
    
    # filter to the same set of genes that supports binning:
    gene_mask = adata.var.index.isin(adata.uns['SCDRS_PARAM']['GENE_STATS'].index)
    print(f"Retaining {np.sum(gene_mask)} / {len(adata.var)} genes")
    adata._inplace_subset_var(gene_mask)
    
    # verbose to print memory usage:
    peak=tracker.stop()
    print(f'peak memory after regression: {peak:.2f} GB') if verbose else None
    return adata if copy else None

def compute_stats(
    adata,
    weight_option,
    ctrl_match_key,
    implicit_cov_corr=False,
    cell_weight=None,
    n_mean_bin=10,
    n_var_bin=10,
    n_length_bin=10,
    n_chunk=20,
    verbose = True,
):
    """
    Compute gene-level and cell-level statstics used for scDRS analysis. `adata`
    should be log-scale. It has two modes. In the normal mode, it computes
    statistics for `adata.X`. In the implicit covariate correction mode, the
    covariate correction has not been performed on `adata.X` but the corresponding
    information is stored in `adata.uns["SCDRS_PARAM"]`. In this case, it computes
    statistics for the covariate-corrected data

        `transformed_X = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed to be log-scale.
    weight_opt : str
        weighting option, same as cli input
    ctrl_match_key : str, default="mean_var"
        Gene-level statistic used for matching control and disease genes;
        should be in `data.uns["SCDRS_PARAM"]["GENE_STATS"]`.
    implicit_cov_corr : bool, default=False
        If True, compute statistics for the implicit corrected data
        `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`. Otherwise, compute
        for the original data `adata.X`.
    cell_weight : array_like, default=None
        Cell weights of length `adata.shape[0]` for cells in `adata`,
        used for computing weighted gene-level statistics.
    n_mean_bin : int, default=10
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=10
        Number of expression-variance bins for matching control genes.
    n_length_bin : int, default=10
        Number of gene-length bins for matching control genes.
    n_chunk : int, default=10
        Number of chunks to split the data into when computing mean and variance
        using _get_mean_var_implicit_cov_corr.
    verbose: bool, default=True

    Returns
    -------
    df_gene : pandas.DataFrame
        Gene-level statistics of shape (n_gene, 7):
        
        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.
        - "length" : gene length.
        - "var_tech" : technical variance in log scale if weight option is vs
        - "ct_mean" : mean expression in original non-log scale. if weight option is vs
        - "ct_var" : variance expression in original non-log scale if weight option is vs
        - "ct_var_tech" : technical variance in original non-log scale if weight option is vs
        - ctrl_match_key : n_mean_bin * n_var_bin mean-variance bins, noted as mean_bin.var_bin the gene belongs to
        
    df_cell : pandas.DataFrame
        Cell-level statistics of shape (n_cell, 2):

        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.
    """

    if implicit_cov_corr:
        assert (
            "SCDRS_PARAM" in adata.uns
        ), "adata.uns['SCDRS_PARAM'] is not found, run `scdrs.pp.preprocess` before calling this function"

    df_gene = pd.DataFrame(
        index=adata.var_names,
        columns=[
            "mean",
            "var",
            "length",
            "var_tech",
            "ct_mean",
            "ct_var",
            "ct_var_tech",
            ctrl_match_key,
        ],
    )
    df_cell = pd.DataFrame(index=adata.obs_names, columns=["mean", "var"])
    
    # Gene-level statistics
    # put the gene length into the dataframe:
    conversion = met_scdrs.util.load_homolog_mapping('mouse', 'human')
    if (df_gene.index.isin(conversion.keys()).sum() > df_gene.index.isin(conversion.values()).sum()):
        # if mouse gene is being used:
        gene_info = met_scdrs.util.load_gencode('mouse')
    else:
        # if human gene is being used:
        gene_info = met_scdrs.util.load_gencode('human')
    
    # get the gene length 
    gene_lenth_dict = dict(zip(gene_info['gene_name'], gene_info['gene_length']))
    df_gene['length'] = df_gene.index.map(gene_lenth_dict)
    
    if not implicit_cov_corr:
        # measure ram usage:
        tracker = met_scdrs.util.MemoryTracker()
        tracker.start()
        
        # Normal mode
        df_gene["mean"], df_gene["var"] = _get_mean_var(
            adata.X, axis=0, weights=cell_weight
        )
        peak = tracker.stop()
        print(f'peak memory during normal mode of _get_mean_var: {peak:.2f} GB') if verbose else None
        
        # Get the mean and var for the non-log-scale size-factor-normalized counts
        # It is highly correlated to the non-size-factor-normalized counts
        tracker.start()
                
        if weight_option == 'vs':
            if sparse.issparse(adata.X):  # sp sparse matrix
                temp_X = adata.X.copy().expm1()  # exp(X)-1 to get ct matrix from logct
            else:
                temp_X = np.expm1(adata.X)  # numpy ndarray
            df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var(
                temp_X, axis=0, weights=cell_weight
            )
            del temp_X
        
        peak = tracker.stop()
        print(f'peak memory during count mean var: {peak:.2f} GB') if verbose else None
    
    else:
        # Implicit covariate correction mode
        df_gene["mean"], df_gene["var"] = _get_mean_var_implicit_cov_corr(
            adata, axis=0, n_chunk=n_chunk, weights=cell_weight
        )
        if weight_option == 'vs':
            df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var_implicit_cov_corr(
                adata, transform_func=np.expm1, axis=0, n_chunk=n_chunk, weights=cell_weight
            )

    tracker.start()
    
    if weight_option == 'vs':    
        # Borrowed from scanpy _highly_variable_genes_seurat_v3
        not_const = df_gene["ct_var"].values > 0
        if (df_gene["ct_mean"].values <= 0).sum() > 0:
            # Exclude genes with negative values (usually small)
            warnings.warn("%d genes with ct_mean<0" % (df_gene["ct_mean"].values < 0).sum())
            not_const = not_const & (df_gene["ct_mean"].values > 0)
        estimat_var = np.zeros(adata.shape[1], dtype=np.float64)
        y = np.log10(df_gene["ct_var"].values[not_const])
        x = np.log10(df_gene["ct_mean"].values[not_const])
        model = loess(x, y, span=0.3, degree=2)
        model.fit()
        estimat_var[not_const] = model.outputs.fitted_values
        df_gene["ct_var_tech"] = 10 ** estimat_var
        # Recipe from Frost Nucleic Acids Research 2020
        df_gene["var_tech"] = df_gene["var"] * df_gene["ct_var_tech"] / df_gene["ct_var"]
        df_gene.loc[df_gene["var_tech"].isna(), "var_tech"] = 0

    peak = tracker.stop()
    print(f'peak memory during tech_var computation: {peak:.2f} GB') if verbose else None
    
    tracker.start()
    
    ###########################################################################################
    ######               rework scDRS control gene matching                              ######
    ###########################################################################################
    bin_param = {'mean': n_mean_bin, 'var': n_var_bin, 'length': n_length_bin}
    bin_counts = {k: bin_param[k] for k in ctrl_match_key.split('_') if k in bin_param}
    
    # if require matching by length, remove the genes that we cnanot find length info
    if 'length' in bin_counts.keys():
        print(f"{df_gene['length'].isna().sum()} genes are removed because no length available")
        df_gene = df_gene.loc[df_gene['length'].notna(), :]
    
    # assign gene to bins:
    df_gene.loc[:, ctrl_match_key] = _assign_gene_to_bins(df_gene, ctrl_match_key, bin_counts)
    
    # print message on memory trace:
    peak = tracker.stop()
    print(f'peak memory during bin computation: {peak:.2f} GB') if verbose else None
    
    tracker.start()
    # Cell-level statistics
    if not implicit_cov_corr:
        # Normal mode
        df_cell["mean"], df_cell["var"] = _get_mean_var(adata.X, axis=1)
    else:
        # Implicit covariate correction mode
        df_cell["mean"], df_cell["var"] = _get_mean_var_implicit_cov_corr(
            adata, axis=1, n_chunk=n_chunk
        )
    
    peak = tracker.stop()
    print(f'peak memory during cell stat computation: {peak:.2f} GB') if verbose else None
    
    return df_gene, df_cell





##############################################################################
######################## Preprocessing Subroutines ###########################
##############################################################################
def _assign_gene_to_bins(df, ctrl_match_key, bin_counts):
    """
    Dynamically use qcut to assign bins to genes:
    
    Parameters
    ----------
    df : pd.DataFrame
        df_gene Dataframe in the preprocess function
    ctrl_match_key : str
        axis to match control genes with, must be separated by '_', e.g.: mean_var
    bin_counts : dictionary
        dictioary where keys are matchable axis, and value is the number of bins to specify
        
    """
    # assertion:
    assert set(bin_counts.keys()) == set(ctrl_match_key.split('_')), "bin_counts keys do not match ctrl_match_key"
    
    features = ctrl_match_key.split('_')
    df = df.copy()
    bin_label_cols = []

    for feature in features:
        new_bin = pd.Series(pd.NA, index=df.index, dtype="Int64")

        # Iterate through unique combinations of prior bins
        grouped = df.groupby(bin_label_cols or [lambda x: 0])  # dummy group if first

        for name, subset in grouped:
            valid = subset[feature].notna()
            if valid.sum() == 0:
                continue
            try:
                binned = pd.qcut(
                    subset.loc[valid, feature],
                    q=bin_counts[feature],
                    labels=False,
                    duplicates='drop'
                ).astype("Int64")
                new_bin.loc[valid.index] = binned
            except ValueError:
                # If not enough values to bin, keep NA
                pass

        col_name = f"{feature}_bin"
        df[col_name] = new_bin
        bin_label_cols.append(col_name)

    # Combine bins into final string label
    df['bin_label'] = df[bin_label_cols].apply(lambda row: '.'.join(row.dropna().astype(str)), axis=1)
    return df['bin_label']
        
def reg_out(mat_Y, mat_X):
    """Regress mat_X out of mat_Y.

    Parameters
    ----------
    mat_Y : np.ndarray
        Response variable of shape (n_sample, n_response).
    mat_X : np.ndarray
        Covariates of shape (n_sample, n_covariates).

    Returns
    -------
    mat_Y_resid : np.ndarray
        Response variable residual of shape (n_sample, n_response).
    """

    if sparse.issparse(mat_X):
        mat_X = mat_X.toarray()
    else:
        mat_X = np.array(mat_X)
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])

    if sparse.issparse(mat_Y):
        mat_Y = mat_Y.toarray()
    else:
        mat_Y = np.array(mat_Y)
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    n_sample = mat_Y.shape[0]
    mat_xtx = np.dot(mat_X.T, mat_X) / n_sample
    mat_xty = np.dot(mat_X.T, mat_Y) / n_sample
    mat_coef = np.linalg.solve(mat_xtx, mat_xty)
    mat_Y_resid = mat_Y - mat_X.dot(mat_coef)

    if mat_Y_resid.shape[1] == 1:
        mat_Y_resid = mat_Y_resid.reshape([-1])
    
    # if the original float is 32, keep the float 32 residuals as well to prevent up-scaling
    if mat_Y.dtype == 'float32':
        return mat_Y_resid.astype(np.float32)
    else:
        return mat_Y_resid

def _reg_out_inplace(mat_Y, mat_X):
    """Regress mat_X out of mat_Y.
    
    Parameters
    ----------
    mat_Y : np.ndarray
        Response variable of shape (n_sample, n_response).
    mat_X : np.ndarray
        Covariates of shape (n_sample, n_covariates).

    Returns
    -------
    mat_Y : np.ndarray
        Response variable residual of shape (n_sample, n_response).
    """

    if sparse.issparse(mat_X):
        mat_X = mat_X.toarray()
    else:
        mat_X = np.asarray(mat_X)
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])
    
    if sparse.issparse(mat_Y):
        mat_Y = mat_Y.toarray()
    else:
        mat_Y = np.asarray(mat_Y)
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])
    
    # ensure float 32 consistency:
    mat_X = mat_X.astype(np.float32, copy = False)
    mat_Y = mat_Y.astype(np.float32, copy = False)
    
    # Compute inverse(X_t * X) * X_t with scaling trick:
    n_sample = mat_Y.shape[0]
    mat_xtx = np.dot(mat_X.T, mat_X) / n_sample
    
    # regress in place:
    for j in tqdm(range(mat_Y.shape[1])):
        y_j = mat_Y[:, j]
        xty_j = np.dot(mat_X.T, y_j) / n_sample
        coef_j = np.linalg.solve(mat_xtx, xty_j)
        mat_Y[:, j] -= np.dot(mat_X, coef_j)
    
    if mat_Y.shape[1] == 1:
        mat_Y = mat_Y.reshape([-1])
    
    # if the original float is 32, keep the float 32 residuals as well to prevent up-scaling
    if mat_Y.dtype != 'float32':
        return mat_Y.astype(np.float32)
    else:
        return mat_Y

def _get_mean_var(sparse_X, axis=0, weights=None):
    """
    Compute mean and var of a sparse / non-sparse matrix.

    Parameters
    ----------
    sparse_X : array_like
        Data matrix (can be dense/sparse).
    axis : {0, 1}, default=0
        Axis along which to compute mean and variance.
    weights : array_like, default=None
        Weights of length `sparse_X.shape[axis]`.

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    v_var : np.ndarray
        Variance vector.

    """

    if weights is None:
        if sparse.issparse(sparse_X):
            # Sparse + unweighted
            v_mean = sparse_X.mean(axis=axis)
            v_mean = np.array(v_mean).reshape([-1])
            v_var = sparse_X.power(2).mean(axis=axis)
            v_var = np.array(v_var).reshape([-1])
            v_var = v_var - v_mean ** 2
        else:
            # Dense + unweighted
            sparse_X = np.array(sparse_X)
            v_mean = np.mean(sparse_X, axis=axis)
            v_var = np.var(sparse_X, axis=axis)

    else:
        weights = np.array(weights)
        if sparse.issparse(sparse_X):
            # Sparse + weighted
            v_mean = _weighted_sparse_average(sparse_X, weights, axis=axis)
            v_var = _weighted_sparse_average(sparse_X.power(2), weights, axis=axis)
            v_var = v_var - v_mean ** 2
        else:
            # Dense + weighted
            sparse_X = np.array(sparse_X)
            v_mean = np.average(sparse_X, axis=axis, weights=weights)
            v_var = np.average(np.square(sparse_X), axis=axis, weights=weights)
            v_var = v_var - v_mean ** 2

    return v_mean, v_var


def _weighted_sparse_average(sparse_X, weights, axis=0):
    """
    Compute weighted mean for a sparse matrix.

    Parameters
    ----------
    sparse_X : array_like
        Sparse data matrix.
    weights : array_like
        Weights of length `sparse_X.shape[axis]`.
    axis : {0, 1}, default=0
        Axis along which to compute mean.

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    """

    if axis == 0:
        v_mean = sparse_X.T.multiply(weights).mean(axis=1)
        v_mean = np.array(v_mean).reshape([-1]) / np.mean(weights)
        return v_mean

    if axis == 1:
        v_mean = sparse_X.multiply(weights).mean(axis=1)
        v_mean = np.array(v_mean).reshape([-1]) / np.mean(weights)
        return v_mean


def _get_mean_var_implicit_cov_corr(
    adata, axis=0, weights=None, transform_func=None, n_chunk=20
):
    """
    Compute mean and variance of sparse matrix of the form

        `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`.

    Computed iteratively over chunks of sparse matrix by converting to
    dense matrix and computing mean and variance of dense matrix.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix (n_obs, n_vars)
    axis : {0, 1}, default=0
        Axis along which to compute mean and variance
    weights : array_like, default=None
        Weights of length `adata.shape[axis]`.
    transform_func : function, default=None
        Function to transform the data before computing mean and variance
    n_chunk : int, default=20
        Number of chunks to split the data into when computing mean and variance
        this will determine the memory usage

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    v_var : np.ndarray
        Variance vector.
    """

    assert axis in [0, 1], "axis must be one of [0, 1]"
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `preprocess` before calling this function"

    n_obs, n_gene = adata.shape

    # COV INFO: cov_mat: (n_cell, n_cov) cov_beta (n_cov, n_gene)
    cell_list = list(adata.obs_names)
    cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
    gene_list = list(adata.var_names)
    cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].loc[cell_list, cov_list].values
    cov_beta = adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
    gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

    if transform_func is None:
        transform_func = lambda x: x

    if weights is not None:
        weights = np.array(weights)

    # mean / var of transform_func(X + cov_mat @ cov_beta + mean)
    if axis == 0:
        # output would have the shape of (n_gene, )
        v_mean = np.zeros(n_gene)
        v_var = np.zeros(n_gene)
        start = 0
        chunk_size = n_gene // n_chunk
        while start < n_gene:
            stop = min(start + chunk_size, n_gene)
            chunk_X = (
                adata.X[:, start:stop]
                + cov_mat @ cov_beta[:, start:stop]
                + gene_mean[start:stop]
            )
            chunk_X = transform_func(chunk_X)
            v_mean[start:stop], v_var[start:stop] = _get_mean_var(
                chunk_X, axis=axis, weights=weights
            )
            start = stop

    elif axis == 1:
        # output would have the shape of (n_obs, )
        v_mean = np.zeros(n_obs)
        v_var = np.zeros(n_obs)
        start = 0
        chunk_size = n_obs // n_chunk
        while start < n_obs:
            stop = min(start + chunk_size, n_obs)
            chunk_X = (
                adata.X[start:stop, :] + cov_mat[start:stop] @ cov_beta + gene_mean
            )
            chunk_X = transform_func(chunk_X)
            v_mean[start:stop], v_var[start:stop] = _get_mean_var(
                chunk_X, axis=axis, weights=weights
            )
            start = stop

    return v_mean, v_var