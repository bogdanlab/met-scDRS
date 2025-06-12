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
import re

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
    transformation : str = None,
    cov_file: str = None,
    ctrl_match_opt: str = 'mean_var',
    weight_opt: str = 'inv_std',
    n_ctrl: int = 1000,
    flag_return_ctrl_raw_score: bool = False,
    flag_return_ctrl_norm_score: bool = True,
    diagnostic = False,
    diagnostic_dir = None,
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
    transformation : str, optional
        "logit" : logit transformation
        "arcsine" : arcsine transformation
        "log_library" : inverse the fraction then library normalize
        if None, no transformation is applied
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
    diagnostic : bool, optional
        if True, plot out the control gene set drawn bins. Default is False.
    diagnostic_dir : str, optional
        directory to store plots for diagnostic. Default is None.
    verbose : bool, optional
        If True, chatty about processing.

    Examples
    --------
    met-scdrs compute_score \
        --h5ad_file <h5ad_file> \
        --preprocess True \
        --preprocess_method inverse \
        --variance_clip 5 \
        --transformation arcsine \
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
        --diagnostic False \
        --diagnostic_dir None \
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
    TRANSFORMATION = transformation
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
    DIAGNOSTIC = diagnostic
    DIAGNOSTIC_DIR = diagnostic_dir
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
    header += "--transformation %s \\\n" % TRANSFORMATION
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
    header += "--diagnostic %s\n" % DIAGNOSTIC
    header += "--diagnostic_dir %s\n" % DIAGNOSTIC_DIR
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
    if CTRL_MATCH_OPT not in ["mean_var_length", "mean_var", "mean", "var", "mean_length", "var_length"]:
        raise ValueError("--ctrl_match_opt mean_var should be one of mean_var_length, mean_var, mean, var, mean_length, var_length")
    if WEIGHT_OPT not in ['inv_std', 'vs']:
        raise ValueError("--weight_opt should be inv_std or vs")
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
            variance_clip = VARIANCE_CLIP,
            transformation = TRANSFORMATION,
            verbose = VERBOSE
        )
        peak = tracker.stop()
        print(f'peak memory during normalization: {peak:.2f} GB') if VERBOSE else None

    # preprocess with covariates
    met_scdrs.preprocess(
        adata,
        cov=df_cov,
        n_mean_bin=10,
        n_var_bin=10,
        n_length_bin = 10,
        copy=False,
        weight_option=WEIGHT_OPT,
        ctrl_match_key=CTRL_MATCH_OPT,
        verbose = VERBOSE)
    
    if DIAGNOSTIC:
        met_scdrs.diagnostic.ctrl_match_bin(
            adata,
            dict_gs = dict_gs,
            ctrl_match_key=CTRL_MATCH_OPT,
            plot_dir = DIAGNOSTIC_DIR)
    print("")
    
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
        print(f'peak memory during score computation: {peak:.2f} GB') if VERBOSE else None
        
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

def compare_score(
    score1 : str,
    score2 : str,
    plot_path : str = None,
):
    """
    score1 : str
        path to first set of disease score file, either a directory, or a score file ends with .score.gz
    score2 : str
        path to second set of disease score file, either a directory, or a score file ends with .score.gz
    plot_path : str
        path to plotting a scatter plot between the two scores, default to None
    """
    # print a header so people know the input:
    header = get_cli_head()
    header += 'Called compare_score \\\n'
    header += "--score1_path %s \\\n" % score1
    header += "--score2_path %s \\\n" % score2
    header += "--plot_path %s \n" %plot_path
    print(header)
    
    # check if the file is directory or file path:
    if os.path.isdir(score1) and os.path.isdir(score2) and os.path.isdir(plot_path):
        # read in the files:
        score1_files = os.listdir(score1)
        score2_files = os.listdir(score2)
        
        # identify the score files:
        score1_files = [score_file for score_file in score1_files if score_file.endswith('.score.gz')]
        score2_files = [score_file for score_file in score2_files if score_file.endswith('.score.gz')]
        common_files = set(score1_files).intersection(score2_files)
        
        # for each file in the common files, make comparisons:
        for f in common_files:
            disease = re.sub('.score.gz', '', f)
            score1_ = os.path.join(score1, f)
            score2_ = os.path.join(score2, f)
            score1_ = pd.read_csv(score1_, sep = '\t', index_col = 0)
            score2_ = pd.read_csv(score2_, sep = '\t', index_col = 0)
            plot_path_ = f"{os.path.join(plot_path, disease)}_comparison.png"
            
            # call the comparison:
            met_scdrs.diagnostic.compare_score(score1_, score2_, plot_path_, add_date = True)
    
    elif os.path.isfile(score1) and os.path.isfile(score2) and os.path.isfile(plot_path):
        # call plotting script:
        score1 = pd.read_csv(score1, sep = '\t', index_col = 0)
        score2 = pd.read_csv(score2, sep = '\t', index_col = 0)
        met_scdrs.diagnostic.compare_score(score1, score2, plot_path, add_date = True)
    
    else:
        print('please check input, if score1, score2, and plot_path can either be all directories, or file paths')

def probe_background(
    score : str,
    plot_path : str,
    sampling : int,
    cell_meta_path : str = None,
    group_column : str = None,
    seed : int = 103
    ):
    """
    investigate the background distribution to see how sampled control scores are distributed
    
    Parameters:
    -----------
    score : str
        path to the score, can either be a directory, or it can be a score file ends with .full_score.gz
    plot_path : str
        path to the plot, can be eitehr a director or a file to plot into
    sampling : int, optional
        if cell_meta and group column provided, sampled in each level, highly recommend to keep it < 30
        if not provided, sampled globally, recommend 100 - 1000 range
        recommendation is made to keep the compute time managable
    cell_group_path : str, optional
        cell meta data path where rows are cells, columns are annotation of cells, default to None
    group_column : str, optional
        group structure of cells to sample background distribution. example: if group_column='cell_type'
        the sampling will be done in each cell type, and it will provide distribution comparison between
        and within cell types, default to None
    seed : int, optional
        numpy sampling seed passed into the sampler, default is 103
    """
    
    # set the cell meta data:
    if cell_meta_path:
        # load in the file:
        cell_meta = pd.read_csv(cell_meta_path, index_col = 0, sep = '\t')
        
        # if the cell group is provided, subset to the pd.Series
        if group_column:
            cell_group = cell_meta.loc[:, group_column]
        else:
            raise AssertionError("please provide cell group column")
    else:
        cell_group = None
    
    # make sure that the score input is a real path or a directory
    assert (os.path.isfile(score) or os.path.isdir(score)), "score needs to be either a file or a directory"
    if os.path.isdir(score) and os.path.isdir(plot_path):
        # grab out the full score files:
        score_files = [file for file in os.listdir(score) if file.endswith('full_score.gz')]
        assert len(score_files) > 0, "there is no full score files, is --flag_return_ctrl_norm_score set to True?"
        
        # for each of these files, do:
        for file in score_files:
            # get disease:
            disease = re.sub('.full_score.gz', '', file)
            # get the score inside:
            score_path = os.path.join(score, file)
            full_score_ = pd.read_csv(score_path, sep = '\t', index_col = 0)
                        
            # get the plot path:
            plot_path_ = os.path.join(plot_path, disease)
            plot_path_ = f"{plot_path_}_calibration_p.png"
            
            # call plot_bg_distribution
            met_scdrs.diagnostic.plot_bg_distribution(full_score_, plot_path_, sampling = sampling, cell_group = cell_group)
    
    elif os.path.isfile(score):
        # directly call the met_scdrs plot bg distribution
        full_score_ = pd.read_csv(score, sep = '\t', index_col = 0)
        met_scdrs.diagnostic.plot_bg_distribution(full_score_, plot_path_, sampling = sampling, cell_group = cell_group)
        
    else:
        print('please check input, score and plot_path can either be all directories, or all file paths')


def quote_from_cyberpunk2077():
    # meant for testing script
    print('To This!')

def main():
    fire.Fire()
