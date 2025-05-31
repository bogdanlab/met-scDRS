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
    cov_file: str = None,
    ctrl_match_opt: str = 'mean_var',
    weight_opt: str = 'inv_std',
    n_ctrl: int = 1000,
    flag_return_ctrl_raw_score: bool = False,
    flag_return_ctrl_norm_score: bool = True
):
    """
    Compute met-scDRS scores for single-cell methylation data.

    Parameters
    ----------
    h5ad_file : str
        Path to the single-cell `.h5ad` file. The `.X` attribute should contain a cell-by-gene 
        methylation matrix.
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

    Examples
    --------
    met-scdrs compute_score \
        --h5ad_file <h5ad_file> \
        --h5ad_species human \
        --gs_file <gs_file> \
        --gs_species human \
        --out_folder <out_folder> \
        --cov_file <cov_file> \
        --ctrl_match_opt mean_var \
        --weight_opt inv_std \
        --n_ctrl 1000 \
        --flag_return_ctrl_raw_score False \
        --flag_return_ctrl_norm_score True
    """
    # record system start time:
    sys_start_time = time.time()

    ###########################################################################################
    ######                                   INPUT PARSEING                              ######
    ###########################################################################################
    H5AD_FILE = h5ad_file
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

    # if the species name doesn't match, convert the species name:
    if H5AD_SPECIES != GS_SPECIES:
        H5AD_SPECIES = met_scdrs.util.convert_species_name(H5AD_SPECIES)
        GS_SPECIES = met_scdrs.util.convert_species_name(GS_SPECIES)
        print(f'h5ad species name converted to {H5AD_SPECIES}')
        print(f'gs_species name converted to {GS_SPECIES}')

    ###########################################################################################
    ######                                    Section Title                              ######
    ###########################################################################################
    






def quote_from_cyberpunk2077():
    # meant for testing script
    print('To This!')

def main():
    fire.Fire()
