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
matplotlib.use('Agg')

def _bin_formatter(gene_ctrol_match_bin):
    """
    Subroutine under ctrl_match_bin that reforats that decondense bins
    
    Parameters
    ----------
    gene_ctrol_match_bin : pd.Series
        essentially: preprocessed_h5ad.uns['SCDRS_PARAM']['GENE_STATS'][ctrl_match_key]
    
    Output
    ------
    decompressed : pd.DataFrame
        if two axis control match (e.g.: mean, and variance), return a heatmap where
        rows are bins for the first axis, and columns are bins for the second axis
        
        if three axis control match (e.g.: mean, variance and gene length), return
        a four column dataframe.
        The first three columns corresponds to the bin membership in each of the three axis. 
        The fourth column (count) represents the number of genes that are in the bin set
        by the three axis.
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
        decompressed = average_gene_in_bin
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
    plots

    
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
            print('do something else')
            plt.figure(figsize=(8, 6))
            g = sns.FacetGrid(grouped_bins, col='bin_axis2', col_wrap=5, height=4)
            g.map_dataframe(
                lambda data, color: sns.heatmap(
                    data.pivot(index='bin_axis0', columns='bin_axis1', values='count').fillna(0),
                    annot=True, fmt='g', cmap='viridis', cbar=False
                )
            )
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
