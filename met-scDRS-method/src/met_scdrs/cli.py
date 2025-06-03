#!/usr/bin/env python

# import packages:
import fire
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import met_scdrs
from typing import Dict, List
import scanpy as sc
import os
import time

### HEADER ########################################################################################
def get_cli_head():
    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* methylation single-cell disease relevance score (met-scDRS)\n"
    MASTHEAD += "* Version %s\n" % met_scdrs.__version__
    MASTHEAD += "* Xinzhe Li\n"
    MASTHEAD += "* Adapted from scDRS written by Martin Jinye Zhang and Kangcheng Hou\n"
    MASTHEAD += "* UPENN\n"
    MASTHEAD += "* MIT License\n"
    MASTHEAD += "******************************************************************************\n"
    return MASTHEAD
    
def compute_score(
    h5ad_file: str,
    h5ad_species: str,
    gs_file: str,
    gs_species: str,
    out_folder: str,
    preprocess: bool = True,
    preprocess_method: str = 'inverse',
    variance_clip: int = 5,
    cov_file: str = None,
    ctrl_match_opt: str = 'mean_var',
    weight_opt: str = 'inv_std',
    n_ctrl: int = 1000,
    flag_return_ctrl_raw_score: bool = False,
    flag_return_ctrl_norm_score: bool = True,
    verbose: bool = True
):
    """
    Compute met-scDRS scores for single-cell methylation data.

    Parameters
    ----------
    h5ad_file : str
        Path to the single-cell `.h5ad` file. The `.X` attribute should contain a cell-by-gene 
        methylation matrix.
    preprocess : bool
        if the fraction matrix in `.h5ad` be preprocessed
    preprocess_method : str
        methodology to preprocess the single cell methylation matrix, supported methods:
        "inverse" : inverse the fraction into 1 - X
    variance_clip : int
        only genes with greater than specified percentile will be retained, default to 5th percentile
    h5ad_species : str
        Species of the cells in the `h5ad_file`. Supports automatic gene name translation 
        between human and mouse. As long as the species matches between `h5ad_file` and `gs_file`, 
        the tool will function properly.
    gs_file : str
        Path to the gene set file (same format as scDRS `.gs` files).
    gs_species : str
        Species of the genes in `gs_file`. Supports automatic gene name translation 
        for either human or mouse.
    out_folder : str
        Directory where output files will be saved.
    cov_file : str, optional
        Path to a design matrix for covariate adjustment using a linear model. Format should match 
        that used in scDRS.
    ctrl_match_opt : str, optional
        Option for control gene matching. Currently supports "mean_var".
    weight_opt : str, optional
        Option for technical noise weighting. Supported values: "vs", "inv_std". Default is "inv_std".
    n_ctrl : int, optional
        Number of control gene sets to sample when generating the null distribution. Default is 1000.
    flag_return_ctrl_raw_score : bool, optional
        If True, the full set of raw control scores will be output. Default is False.
    flag_return_ctrl_norm_score : bool, optional
        If True, normalized control scores will be output. Default is True.
    verbose : bool, optional
        If True, chatty about processing

    Examples
    --------
    met-scdrs compute_score \
        --h5ad_file <h5ad_file> \
        --preprocess True \
        --preprocess_method inverse \
        --variance_clip 5 \
        --h5ad_species human \
        --gs_file <gs_file> \
        --gs_species human \
        --out_folder <out_folder> \
        --cov_file <cov_file> \
        --ctrl_match_opt mean_var \
        --weight_opt inv_std \
        --n_ctrl 1000 \
        --flag_return_ctrl_raw_score False \
        --flag_return_ctrl_norm_score True \
        --verbose True
    """
    # record system start time:
    sys_start_time = time.time()
    
    ###########################################################################################
    ######                                   INPUT PARSEING                              ######
    ###########################################################################################
    H5AD_FILE = h5ad_file
    PREPROCESS = preprocess
    PREPROCESS_METHOD = preprocess_method
    VARIANCE_CLIP = variance_clip
    H5AD_SPECIES = h5ad_species
    COV_FILE = cov_file
    GS_FILE = gs_file
    GS_SPECIES = gs_species
    OUT_FOLDER = out_folder
    CTRL_MATCH_OPT = ctrl_match_opt
    WEIGHT_OPT = weight_opt
    N_CTRL = n_ctrl
    FLAG_RETURN_CTRL_RAW_SCORE = flag_return_ctrl_raw_score
    FLAG_RETURN_CTRL_NORM_SCORE = flag_return_ctrl_norm_score
    VERBOSE = verbose
    
    # if the species name doesn't match, convert the species name:
    if H5AD_SPECIES != GS_SPECIES:
        H5AD_SPECIES = met_scdrs.util.convert_species_name(H5AD_SPECIES)
        GS_SPECIES = met_scdrs.util.convert_species_name(GS_SPECIES)
        print(f'h5ad species name converted to {H5AD_SPECIES}')
        print(f'gs_species name converted to {GS_SPECIES}')
    
    # print out the header:
    header = get_cli_head()
    header += "Call: met-scdrs compute-score \\\n"
    header += "--h5ad_file %s \\\n" % H5AD_FILE
    header += "--preprocess %s \\\n" % PREPROCESS
    header += "--preprocess_method %s \\\n" % PREPROCESS_METHOD
    header += "--variance_clip %s \\\n" % VARIANCE_CLIP
    header += "--h5ad_species %s \\\n" % H5AD_SPECIES
    header += "--cov_file %s \\\n" % COV_FILE
    header += "--gs_file %s \\\n" % GS_FILE
    header += "--gs_species %s \\\n" % GS_SPECIES
    header += "--ctrl_match_opt %s \\\n" % CTRL_MATCH_OPT
    header += "--weight_opt %s \\\n" % WEIGHT_OPT
    header += "--n_ctrl %d \\\n" % N_CTRL
    header += "--flag_return_ctrl_raw_score %s \\\n" % FLAG_RETURN_CTRL_RAW_SCORE
    header += "--flag_return_ctrl_norm_score %s \\\n" % FLAG_RETURN_CTRL_NORM_SCORE
    header += "--out_folder %s\n" % OUT_FOLDER
    header += "--verbose %s \n" % VERBOSE
    print(header)
    
    # check options:
    print('\n\n CHECKING INPUT \n\n')
    if H5AD_SPECIES != GS_SPECIES:
        if H5AD_SPECIES not in ["mmusculus", "hsapiens"]:
            raise ValueError(
                "--h5ad-species needs to be one of [mmusculus, hsapiens] "
                "unless --h5ad-species==--gs-species"
            )
        if GS_SPECIES not in ["mmusculus", "hsapiens"]:
            raise ValueError(
                "--gs-species needs to be one of [mmusculus, hsapiens] "
                "unless --h5ad-species==--gs-species"
            )
    # matching control genes should be mean_var:
    if CTRL_MATCH_OPT not in ["mean_var"]:
        raise ValueError("--ctrl_match_opt mean_var should be mean_var")
    if WEIGHT_OPT not in ['inv_std']:
        raise ValueError("--weight_opt should be inv_std")
    # also check folder:
    if not met_scdrs.util.check_folder(OUT_FOLDER):
        raise ValueError("--out_folder does not exist")
    
    ###########################################################################################
    ######                                    DATA LOADING                               ######
    ###########################################################################################
    print('LOADING DATA')
    
    # load h5ad data:
    adata = met_scdrs.util.load_h5ad(H5AD_FILE)
    print(
        "--h5ad-file loaded: n_cell=%d, n_gene=%d (sys_time=%0.1fs)"
        % (adata.shape[0], adata.shape[1], time.time() - sys_start_time)
    )
    print("First 3 cells: %s" % (str(list(adata.obs_names[:3]))))
    print("First 5 genes: %s" % (str(list(adata.var_names[:5]))))
    
    # load in covariate file:
    # Load .cov file
    if COV_FILE is not None:
        df_cov = pd.read_csv(COV_FILE, sep="\t", index_col=0)
        df_cov.index = [str(x) for x in df_cov.index]
        print(
            "--cov-file loaded: covariates=%s (sys_time=%0.1fs)"
            % (str(list(df_cov.columns)), time.time() - sys_start_time)
        )
        print(
            "n_cell=%d (%d in .h5ad)"
            % (df_cov.shape[0], len(set(df_cov.index) & set(adata.obs_names)))
        )
        print("First 3 cells: %s" % (str(list(df_cov.index[:3]))))
        for col in df_cov.columns:
            print(
                "First 5 values for '%s': %s" % (col, str(list(df_cov[col].values[:5])))
            )
    else:
        df_cov = None
    
    # LOAD gs file:
    dict_gs = met_scdrs.util.load_gs(
        GS_FILE,
        src_species=GS_SPECIES,
        dst_species=H5AD_SPECIES,
        to_intersect=adata.var_names,
    )
    print(
        "--gs-file loaded: n_trait=%d (sys_time=%0.1fs)"
        % (len(dict_gs), time.time() - sys_start_time)
    )
    print("Print info for first 3 traits:")
    for gs in list(dict_gs)[:3]:
        print(
            "First 3 elements for '%s': %s, %s"
            % (gs, str(dict_gs[gs][0][:3]), str(dict_gs[gs][1][:3]))
        )
    print("")
    
    # print a message on total time spent
    print(f'loading completed, time elapsed: {time.time() - sys_start_time:.3f} seconds')
    
    ###########################################################################################
    ######                                    PREPROCESS                                 ######
    ###########################################################################################
    tracker = met_scdrs.util.MemoryTracker()
    if PREPROCESS:
        tracker.start()
        adata = met_scdrs.normalize(
            h5ad_obj = adata,
            method = PREPROCESS_METHOD,
            variance_clip = VARIANCE_CLIP
        )
        peak = tracker.stop()
        if VERBOSE:
            print(f'peak memory during normalization: {peak:.2f} GB')
        
    # preprocess with covariates
    tracker.start()
    met_scdrs.preprocess(adata, cov=df_cov, n_mean_bin=20, n_var_bin=20, copy=False)
    peak = tracker.stop()
    
    # output for testing:
    met_scdrs.util.write_adata_to_csv(adata, '/u/scratch/l/lixinzhe/revision_scratch/batch_regress/scdrs_regress_inplace_fx.csv', subset = True)
    
    if VERBOSE:
        print(f'peak memory during preprocessing: {peak:.2f} GB')
    print("")
    return
    
    ###########################################################################################
    ######                                    Compute score                              ######
    ###########################################################################################
    # Compute score
    print("Computing scDRS score:")
    for trait in dict_gs:
        gene_list, gene_weights = dict_gs[trait]
        if len(gene_list) < 10:
            print(
                "trait=%s: skipped due to small size (n_gene=%d, sys_time=%0.1fs)"
                % (trait, len(gene_list), time.time() - sys_start_time)
            )
            continue
        
        # compute the score (both disease and control genes)
        tracker.start()
        df_res = met_scdrs.score_cell(
            adata,
            gene_list,
            gene_weight=gene_weights,
            ctrl_match_key=CTRL_MATCH_OPT,
            n_ctrl=N_CTRL,
            weight_opt=WEIGHT_OPT,
            return_ctrl_raw_score=FLAG_RETURN_CTRL_RAW_SCORE,
            return_ctrl_norm_score=FLAG_RETURN_CTRL_NORM_SCORE,
            verbose=VERBOSE,
        )
        peak = tracker.stop()
        if VERBOSE:
            print(f'peak memory during score computation: {peak:.2f} GB')
        
        df_res.iloc[:, 0:6].to_csv(
            os.path.join(OUT_FOLDER, "%s.score.gz" % trait),
            sep="\t",
            index=True,
            compression="gzip",
        )
        
        if FLAG_RETURN_CTRL_RAW_SCORE | FLAG_RETURN_CTRL_NORM_SCORE:
            df_res.to_csv(
                os.path.join(OUT_FOLDER, "%s.full_score.gz" % trait),
                sep="\t",
                index=True,
                compression="gzip",
            )
        
        # compute the fdr correction:
        v_fdr = multipletests(df_res["pval"].values, method="fdr_bh")[1]
        n_rej_01 = (v_fdr < 0.1).sum()
        n_rej_02 = (v_fdr < 0.2).sum()
        print(
            "Trait=%s, n_gene=%d: %d/%d FDR<0.1 cells, %d/%d FDR<0.2 cells (sys_time=%0.1fs)"
            % (
                trait,
                len(gene_list),
                n_rej_01,
                df_res.shape[0],
                n_rej_02,
                df_res.shape[0],
                time.time() - sys_start_time,
            )
        )
    return


def quote_from_cyberpunk2077():
    # meant for testing script
    print('To This!')

def main():
    fire.Fire()
