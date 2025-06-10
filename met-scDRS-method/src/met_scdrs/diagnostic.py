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

# define helper:
def _compute_pair(i, j, score_array):
    w_dist = sp.stats.wasserstein_distance(score_array[i], score_array[j])
    return i, j, w_dist

def _plot_calibration(observed, theoretical, plot_path, **kwargs):
    # sort the observed and theoretical pvals:
    sort_observe = np.sort(observed)
    sort_theoretical = np.sort(theoretical)
    
    # make figure;
    plt.figure(figsize=(6, 6))
    plt.plot(-np.log10(sort_theoretical), -np.log10(sort_observe), 'o', markersize=2)
    plt.plot([0, -np.log10(min(sort_theoretical[0], sort_observe[0]))],
            [0, -np.log10(min(sort_theoretical[0], sort_observe[0]))],
            'r--', label='Expected = Observed')
    
    plt.xlabel('Expected -log10(p) (sampled)')
    plt.ylabel('Observed -log10(p)')
    plt.title(f'QQ Plot of P-values (Sampled Null) {kwargs.get("title")}')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_bg_distribution(score, plot_path, sampling = 100, cell_group = None, seed = 103):
    """
    visualize background distribution as a distribution 
    
    Parameters
    ----------
    score : pd.DataFrame 
        a dataframe of full score, if both raw score and control score is present, we plot out both distributions.
    plot_path : str
        a file path where the plot will be outputted to.
    sampling : int
        number of cells to sample from the whole dataframe
    cell_group : pd.Series
        a n_cell by 1 pandas series that documents the cell group identity, default is None
        if the cell_group is provided, sampling is done within each group
    seed : int
        the sampling random seed
    
    Outputs
    -------
    graph1 : a plot of flattened pairwise wesserstein distance distribution between sampled cells
    graph2 : if the cell group is provided, visualize the distribution within group and between group
    graph3 : if sampling is < 10, also plot out the actual contrl norm and raw score for each cell  
    """
    # get the columns that are control row score and control norma score
    control_norm_score_columns = score.columns.str.contains('ctrl_norm_score')
    control_raw_score_columns = score.columns.str.contains('ctrl_raw_score')
    np.random.seed(seed)
    
    # first we will do sampling based on if the cell group is provided or not:
    if cell_group is not None:
        # assert the cell group is the same order as the score
        assert (cell_group.index == score.index).all()
        sampled_cells = []
        # sample from every group
        for group in cell_group.unique():
            group_index = cell_group.index[cell_group == group]
            
            # explicitly check if the group have enough to be sampled:
            if len(group_index) < sampling:
                print(f'{group} has less samples than specified sampling instance, sampling all cells in the group instead')
                sampled_in_group = np.random.choice(group_index, size = len(group_index), replace = False)
            else:
                sampled_in_group = np.random.choice(group_index, size = sampling, replace = False)
            sampled_cells.append(sampled_in_group)
        
        # flatten the sampledcells:
        sampled_cells = np.concatenate(sampled_cells)
    else:
        sampled_cells = np.random.choice(score.index, size = sampling, replace = False)
    
    ###########################################################################################
    ######                            H0: same in shape                                  ######
    ###########################################################################################
    # normalize to zscores:
    if control_raw_score_columns.sum() > 0 and control_norm_score_columns.sum() > 0:
        control_raw_score = score.loc[:, control_raw_score_columns]
        normalized_raw_score = score.loc[:, control_norm_score_columns]
        zscored = scipy.stats.zscore(control_raw_score, axis = 1)
        
        # set the mean cell as a reference for both the raw and raw score
        pvals = []
        norm_pvals = []
        reference = zscored.iloc[0, :]
        norm_ref = normalized_raw_score.iloc[0, :]
        wasserstein_distance = []
        
        for i in tqdm(range(1, len(zscored))):
            stat, pval = scipy.stats.ks_2samp(reference, zscored.iloc[i, :])
            pvals.append(pval)
            stat, pval = scipy.stats.ks_2samp(norm_ref, zscored.iloc[i, :])
            norm_pvals.append(pval)
            wasserstein_distance.append(sp.stats.wasserstein_distance(control_raw_score.iloc[1, :], control_raw_score.iloc[i, :]))
        
        # sort the observed p value
        observed_pval = np.array(pvals)
        sorted_norm_p = np.sort(np.array(norm_pvals))
        sorted_pvals = np.sort(observed_pval)
        
        # draw and sort p value from U(0,1)
        theoretical_p = np.sort(np.random.uniform(0, 1, size = len(sorted_pvals)))
        
        # QQ plot (log scale)
        if plot_path:
            # plot qq plot:
            plot_path_raw = re.sub('.png', '_raw.png', plot_path)
            plot_path_norm = re.sub('.png', '_norm.png', plot_path)
            _plot_calibration(observed_pval, theoretical_p, plot_path_raw, title = 'raw')
            _plot_calibration(norm_pvals, theoretical_p, plot_path_norm, title = 'normalized')
            
            # plot density:
            plot_path_density = re.sub('.png', '_wasserstein_distance_density.png', plot_path)
            sns.kdeplot(wasserstein_distance, fill = True)
            plt.ylabel("density")
            plt.xlabel("control raw score wasserstein distance")
            plt.savefig(plot_path_density, dpi=300, bbox_inches='tight')
            plt.close()
        
        ###########################################################################################
        ######                       plot sampled control score                              ######
        ###########################################################################################
        # also directly visualize 10 of the distributions:
        n_plots = len(sampled_cells)
        if n_plots > 25:
            plots_per_row = 15
        else:
            plots_per_row = 5
        nrows = (n_plots + plots_per_row - 1) // plots_per_row
        
        # initiate figure:
        fig, axes = plt.subplots(nrows, plots_per_row, figsize=(4 * plots_per_row, 3 * nrows))
        axes = axes.flatten()
        
        if plot_path:
            plot_path_raw = re.sub('.png', '_raw_sampling.png', plot_path)
            
            # for sampled cells, plot the density of the control distribution in the sampling:
            for plot_i, sampled_cell in enumerate(sampled_cells):
                cell_type = cell_group.loc[sampled_cell, ]
                values = control_raw_score.loc[sampled_cell, :].values
                sns.kdeplot(values, fill = True, ax = axes[plot_i])
                axes[plot_i].set_ylabel("density")
                axes[plot_i].set_xlabel("control raw score")
                axes[plot_i].set_title(f"{sampled_cell}\n{cell_type}")
            
            for ax in axes[n_plots:]:
                ax.set_visible(False)
            plt.subplots_adjust(hspace=0.8, wspace=0.3)
            plt.savefig(plot_path_raw, dpi=300, bbox_inches='tight')
            plt.close()
        
        # do the same for the normalized control scores:
        fig, axes = plt.subplots(nrows, plots_per_row, figsize=(4 * plots_per_row, 3 * nrows))
        axes = axes.flatten()
        
        if plot_path:
            plot_path_norm = re.sub('.png', '_norm_sampling.png', plot_path)
            
            # for sampled cells, plot the density of the control distribution:
            for plot_i, sampled_cell in enumerate(sampled_cells):
                values = normalized_raw_score.loc[sampled_cell, :].values
                sns.kdeplot(values, fill = True, ax = axes[plot_i])
                axes[plot_i].set_ylabel("density")
                axes[plot_i].set_xlabel("control raw score")
                axes[plot_i].set_title(f"{sampled_cell}\n{cell_type}")
            
            for ax in axes[n_plots:]:
                ax.set_visible(False)
            plt.subplots_adjust(hspace=0.8, wspace=0.3)
            plt.savefig(plot_path_norm, dpi=300, bbox_inches='tight')
            plt.close()
    else:
        print('make sure --flag_return_ctrl_raw_score and --flag_return_ctrl_norm_score is true')
