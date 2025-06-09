from pathlib import Path
import numpy as np
import scipy as sp
import scipy.sparse
import pandas as pd
import numbers
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import os
import anndata
import matplotlib.transforms as mtrans
import matplotlib
import matplotlib.pyplot as plt
from typing import List, Dict
from anndata import read_h5ad
import fnmatch
import matplotlib.patches as patches
import psutil
import threading
import time
import re
import itertools
from tqdm import tqdm
from joblib import Parallel, delayed

matplotlib.use('Agg')

from datetime import date
today = date.today()

def _bin_formatter(gene_ctrol_match_bin):
    """
    Subroutine under ctrl_match_bin that reforats that decondense bins
    
    Parameters
    ----------
    gene_ctrol_match_bin : pd.Series
        essentially: preprocessed_h5ad.uns['SCDRS_PARAM']['GENE_STATS'][ctrl_match_key]
    
    Output
    ------
    axis_num : int
        number of axis in matching scheme detected
    average_gene_in_bin : float
        the average number of genes in a bin
    decompressed : pd.DataFrame
        if two axis control match (e.g.: mean, and variance), return a wide matrix where
        rows are bins for the first axis(mean), and columns are bins for the second axis (variance)
        
        if three axis control match (e.g.: mean, variance and gene length), return
        dictionary where the keys reporesents the bins in the third axis (gene length)
        the value of the dictionary is the same as two axis control match. i.e. a slice
    """
    # split the gene control bin by binning axis:
    splitted_bins = gene_ctrol_match_bin.str.split('.', expand = True)
    split_bins = [f'bin_axis{num}' for num in range(splitted_bins.shape[1])]
    splitted_bins.columns = split_bins
    
    # group them by splitting axis:
    grouped = splitted_bins.groupby(split_bins)
    
    # judge cases:
    if (len(split_bins) == 1):
        print('not supported yet with single axis binning')
        return
    
    elif len(split_bins) == 2:
        axis_num = len(split_bins)
        print('detected two axis style binning')
        grouped_bins = grouped.size()
        bin_matrix = grouped_bins.unstack(fill_value = 0)
        
        average_gene_in_bin = bin_matrix.mean().mean()
        print(f'average genes in bin {average_gene_in_bin}')
        
        # sort the bins so they look nicer:
        bin_matrix.index = bin_matrix.index.astype(int)
        bin_matrix.columns = bin_matrix.columns.astype(int)
        bin_matrix = bin_matrix.sort_index(axis=0).sort_index(axis=1)
        decompressed = bin_matrix
        return axis_num, average_gene_in_bin, decompressed
    
    elif len(split_bins) == 3:
        axis_num = len(split_bins)
        print('detected three axis style binning')
        grouped_bins = grouped.size().reset_index(name = 'count')
        average_gene_in_bin = grouped_bins['count'].sum() / len(grouped_bins)
        print(f'average genes in bin {average_gene_in_bin}')
        
        # initiate:
        decompressed = {}
        
        # for each of the level in the last axis, create binary matrix:
        for last_axis_bin in grouped_bins['bin_axis2']:
            decompressed[f'length_bin_{last_axis_bin}'] = grouped_bins[grouped_bins.bin_axis2 == last_axis_bin].pivot_table(
                columns = 'bin_axis1',
                index = 'bin_axis0',
                values = 'count',
                fill_value = 0
            )
        return axis_num, average_gene_in_bin, decompressed
    
    else:
        print('detected > three axis style binning')
        print('not supported')
        return
         

def ctrl_match_bin(preprocessed_h5ad, dict_gs, ctrl_match_key, plot_dir = None):
    """
    Generate matrix and plot heatmap on disease and control gene in bin
    
    Parameters
    ----------
    preprocessed_h5ad : anndata.AnnData
        preprocessed by met_scdrs.preprocess.preprocess function
    dict_gs : dict
        gene score file loaded by load_gs
    ctrl_match_key : str
        method used to match control genes, same as --ctrl_match_opt input
    plot_dir : str
        output of heatmap showing number of gene in bin. Default to None
    
    Output
    ------
    plots : png
        fraction plots at plot_dir if plot_dir is not None
        the fraction of disease genes and background genes in each bins are visualized
        plot is used to understand if there is over-sampling at any specific bins
    inline messages : str
        if the user did not specify a plot dir, inline messages will be printed
        for every traits in the dict_gs, print the average number of disease genes in bin
        plot_dir is highly recommended for diagnostic purposes
    
    """
    # check if the matching key is in the preprocessed h5ad:
    gene_stats = preprocessed_h5ad.uns['SCDRS_PARAM']['GENE_STATS']
    assert (ctrl_match_key in gene_stats.columns), 'ctrl_match_key not detected'
    
    # obtain the matrix for how the background gene is distributed:
    control_bin = gene_stats[ctrl_match_key]
    axis_num, average_gene_in_bin, bin_decompressed = _bin_formatter(control_bin)
    
    # for each of the gene set, understand where the disease genes are distributed:
    for disease, geneset in dict_gs.items():
        disease_gene_bin = control_bin.loc[control_bin.index.isin(geneset[0]), ]
        axis_num, average_gene_in_bin, disease_bin_decompressed = _bin_formatter(disease_gene_bin)
        
        if axis_num == 2:
            # refill the bins as 0 and inflate it back to original size.
            disease_bin_decompressed = disease_bin_decompressed.reindex(index=range(bin_decompressed.shape[0]), columns=range(bin_decompressed.shape[1]), fill_value=0)
            
            # get the fraction between the disease gene and the rest of the genes:
            fg_bg_ratio = disease_bin_decompressed.divide(bin_decompressed) * 100
            
            # if plot is engaged:
            if plot_dir:
                # grab out the plot figure path:
                plot_path = os.path.join(plot_dir, f'{disease}_foreground_background_percentage_across_bin.png')
                plt.figure(figsize=(8, 6))
                
                # make the heatmap itself:
                g = sns.heatmap(fg_bg_ratio, annot=True, fmt='.1f', annot_kws = {'size': 8}, cmap='Reds', cbar=True, vmin= 0, vmax = 100, linewidths = 0.5, linecolor = 'lightgray', cbar_kws={'label': '% fg/bg'})
                plt.ylabel("Mean bins")
                plt.xlabel("Variance bins")
                plt.tight_layout()
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                plt.close()
            else:
                print('')
                print(f'processing diseases {disease}')
                print(f'{axis_num} way binning detected, average disease gene in bin: {average_gene_in_bin}')
                print(f'recommend setting plot_dir to obtain getting fraction plot if diagnostic is important')
                print('')
                
        if axis_num == 3:
            n_plots = len(disease_bin_decompressed)
            plots_per_row = 4
            nrows = (n_plots + plots_per_row - 1) // plots_per_row
            
            # initiate figure:
            fig, axes = plt.subplots(nrows, plots_per_row, figsize=(4 * plots_per_row, 3 * nrows))
            axes = axes.flatten()
            
            for plot_i, (gene_length_bin, wide_matrix) in enumerate(disease_bin_decompressed.items()):
                # get the background slice:
                background_slice = bin_decompressed[gene_length_bin]
                background_slice.index = background_slice.index.astype(int)
                background_slice.columns = background_slice.columns.astype(int)
                
                # compute the fraction gene between foreground and background:
                wide_matrix.index = wide_matrix.index.astype(int)
                wide_matrix.columns = wide_matrix.columns.astype(int)
                wide_matrix = wide_matrix.reindex(index = range(background_slice.shape[0]), columns = range(background_slice.shape[1]), fill_value=0)
                fg_bg_ratio = wide_matrix.divide(background_slice) * 100
                
                # make the heatmap
                sns.heatmap(fg_bg_ratio, ax = axes[plot_i], annot=True, fmt='.0f', annot_kws = {'size': 8}, cmap='Reds', cbar=True, vmin= 0, vmax = 100, linewidths = 0.5, linecolor = 'lightgray', cbar_kws={'label': '% fg/bg'})
                axes[plot_i].set_ylabel("Mean bins")
                axes[plot_i].set_xlabel("Var bins")
                axes[plot_i].set_title(f'{gene_length_bin}')
            
            # Delete unused axes if any
            for j in range(n_plots, len(axes)):
                fig.delaxes(axes[j])
            
            if plot_dir:
                plot_path = os.path.join(plot_dir, f'{disease}_foreground_background_percentage_across_length_bins.png')
                plt.subplots_adjust(hspace=0.4, wspace=0.3)
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                plt.close()
            else:
                print('')
                print(f'processing diseases {disease}')
                print(f'{axis_num} way binning detected, average disease gene in bin: {average_gene_in_bin}')
                print(f'recommend setting plot_dir to obtain getting fraction plot if diagnostic is important')
                print('')

def compare_score(score1, score2, plot_path, add_date = True):
    """
    Compare between two met_scdrs disease scores
    
    Parameters
    ----------
    score1 : str
        path to first set of disease score file, should end with .score.gz
    score2 : str
        path to second set of disease score file, should end with .score.gz
    plot_path : str
        path to plot scatter plot, default to None
    add_date : bool
        if date should be added after the directory, default to True
    """
    # assertions:
    assert len(score1) == len(score2), "score 1 and score 2 is not the same length"
    assert score1.index.isin(score2.index).all(), "some cells in score 1 not in score 2"
    score1 = score1.loc[score2.index, :]
    
    # compare spearman correlation between normalized scores:
    cor_df = pd.DataFrame({'pval_score1' : score1.pval, 'pval_score2' : score2.pval})
    cor_coef = cor_df.corr(method = 'spearman').iloc[0, 1]
    print(f'spearman correlation between pvalues = {cor_coef:.2f}')
    
    # compare spearman correlation between normalized scores:
    cor_df = pd.DataFrame({'normalized_score1' : score1.norm_score, 'normalized_score2' : score2.norm_score})
    cor_coef = cor_df.corr(method = 'spearman').iloc[0, 1]
    print(f'spearman correlation between scores = {cor_coef:.2f}')
    
    if plot_path:
        # if we would like to add a date to the basename:
        if add_date:
            plot_dir = os.path.dirname(plot_path)
            basename = os.path.basename(plot_path)
            plot_path = os.path.join(plot_dir, f"{today}_{basename}")
        
        # create a scatter plot:
        plt.figure(figsize=(4, 4))
        sns.scatterplot(data = cor_df, x = 'normalized_score1', y = 'normalized_score2')
        abline(cor_df.normalized_score1, cor_df.normalized_score2)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

# define helper:
def _compute_pair(i, j, score_array):
    w_dist = sp.stats.wasserstein_distance(score_array[i], score_array[j])
    return i, j, w_dist

def plot_bg_distribution(score, plot_path):
    """
    visualize background distribution as a distribution 
    
    Parameters
    ----------
    score : pd.DataFrame 
        a dataframe of full score, if both raw score and control score is present, we plot out both distributions.
    plot_path : str
        a file path where the plot will be outputted to. 
    """
    # get the columns that are control row score and control norma score
    control_norm_score_columns = score.columns.str.contains('ctrl_norm_score')
    control_raw_score_columns = score.columns.str.contains('ctrl_raw_score')
    
    # if the control norm score is present:
    if control_norm_score_columns.sum() > 0:
        # build a nested for loop for all cell pairwise distance:
        looper = list(itertools.combinations(range(len(score)), 2))
        looper_length = len(looper)
        
        # optimized for memory:
        score_array = score.loc[:, control_norm_score_columns].values
        results = Parallel(n_jobs = -1, require = 'sharedmem')(delayed(_compute_pair)(arg1, arg2, score_array) for arg1, arg2 in tqdm(looper, total = looper_length))
        breakpoint()
    
def abline(x, y, **kwargs):
    """
    add a line with X=Y
    x : the vector used as x for scatterplot
    y : the vector used as y for scatterplot
    """
    limit_min = min(np.min(x), np.min(y))
    limit_max = max(np.max(x), np.max(y))
    lims = [limit_min, limit_max]
    plt.plot(lims, lims, '--', color = 'black', **kwargs)
