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


###########################################################################################
######                              met-scdrs utilities                              ######
###########################################################################################
def check_folder(folder_path):
    folder = Path(folder_path)
    if folder.exists() and folder.is_dir():
        return True
    else:
        return False

def get_memory():
    return psutil.Process(os.getpid()).memory_info().rss / 1e6

def write_adata_to_csv(h5ad_obj, csv_out, subset = True):
    # Convert adata.X to a dense matrix if it's sparse
    if scipy.sparse.issparse(h5ad_obj.X):
        matrix = h5ad_obj.X.toarray()
    else:
        matrix = h5ad_obj.X
    
    # Create DataFrame with proper labels
    df = pd.DataFrame(
        matrix,
        index=h5ad_obj.obs_names,   # cell IDs as row index
        columns=h5ad_obj.var_names  # gene names as column names
    )
    
    if subset and len(df) > 1000:
        df = df.iloc[:1000, ]
    
    df.to_csv(csv_out)
    return df
    
class MemoryTracker:
    def __init__(self, interval = 0.01):
        self.interval = interval
        self._peak = 0
        self._stop = False
    
    def _track(self):
        process = psutil.Process(os.getpid())
        while not self._stop:
            current = process.memory_info().rss / 1e6 / 1e3
            self._peak = max(self._peak, current)
            time.sleep(self.interval)
    def start(self):
        self._stop = False
        self._peak = 0
        self._thread = threading.Thread(target = self._track)
        self._thread.start()
    def stop(self):
        self._stop = True
        self._thread.join()
        return self._peak



###########################################################################################
######                                  scdrs utilities                              ######
###########################################################################################
def load_h5ad(h5ad_file : str):
    """
    Load in h5ad file and check if there is NA in the data:

    Parameters
    ----------
    h5ad_file : str
        Path to the single-cell `.h5ad` file. The `.X` attribute should contain a cell-by-gene 
        methylation matrix.
    """
    adata = read_h5ad(h5ad_file)
    
    # check inputs (1) no NaN in adata.X
    if np.isnan(adata.X.sum()):
        raise ValueError(
            "h5ad expression matrix should not contain NaN. Please impute them beforehand."
        )
    return adata

def convert_species_name(species):
    if species in ["Mouse", "mouse", "Mus_musculus", "mus_musculus", "mmusculus"]:
        return "mmusculus"
    if species in ["Human", "human", "Homo_sapiens", "homo_sapiens", "hsapiens"]:
        return "hsapiens"
    raise ValueError("species name '%s' is not supported" % species)


def load_homolog_mapping(src_species: str, dst_species: str) -> dict:
    """Load gene homologs between mouse and human.

    Parameters
    ----------
    src_species : str
        Source species. One of 'mmusculus', 'mouse', 'hsapiens', or 'human'.
    dst_species : str
        Destination species. One of 'mmusculus', 'mouse', 'hsapiens', or 'human'.
        Cannot be the same as `src_species`.

    Returns
    -------
    dic_map : dict
        Dictionary of gene homologs (gene symbol).
    """

    src_species = convert_species_name(src_species)
    dst_species = convert_species_name(dst_species)

    assert src_species != dst_species, "src and dst cannot be the same"

    df_hom = pd.read_csv(
        os.path.join(os.path.dirname(__file__), "data/mouse_human_homologs.txt"),
        sep="\t",
    )
    if (src_species == "hsapiens") & (dst_species == "mmusculus"):
        dic_map = {
            x: y for x, y in zip(df_hom["HUMAN_GENE_SYM"], df_hom["MOUSE_GENE_SYM"])
        }
    elif (src_species == "mmusculus") & (dst_species == "hsapiens"):
        dic_map = {
            x: y for x, y in zip(df_hom["MOUSE_GENE_SYM"], df_hom["HUMAN_GENE_SYM"])
        }
    else:
        raise NotImplementedError(
            f"gene conversion from {src_species} to {dst_species} is not supported"
        )

    return dic_map

def load_gencode(species):
    """
    Load gene info from gencode
    
    Parameters
    ----------
    species : str
        One of 'mmusculus', 'mouse', 'hsapiens', or 'human'
    """
    species_name = convert_species_name(species)
    
    if species_name == 'mmusculus':
        gene_info = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "data/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"),
            sep = "\t",
            comment = '#',
            header = None,
            names = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        )
        
        # Keep only gene entries
        gene_info = gene_info[gene_info['feature'] == 'gene']
        gene_info['gene_name'] = gene_info['attribute'].str.extract('gene_name "([^"]+)"')
        gene_info['gene_length'] = gene_info['end'] - gene_info['start'] + 1
    
    else:
        gene_info = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "data/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"),
            sep = "\t",
            comment = '#',
            header = None,
            names = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        )
        
        # Keep only gene entries
        gene_info = gene_info[gene_info['feature'] == 'gene']
        gene_info['gene_name'] = gene_info['attribute'].str.extract('gene_name "([^"]+)"')
        gene_info['gene_length'] = gene_info['end'] - gene_info['start'] + 1
    return gene_info

def load_gs(
    gs_path: str,
    src_species: str = None,
    dst_species: str = None,
    to_intersect: List[str] = None,
) -> dict:
    """Load the gene set file (.gs file).

    Parameters
    ----------
    gs_path : str
        Path to the gene set file with two columns 'TRAIT' and 'GENESET', separated by tab.
        'TRAIT' column contain trait names. 'GENESET' column contain gene names (matching
        expression matrix) and gene weights (for weighted gene set). For unweighted gene set,
        the 'GENESET' column contains only gene names separated by comma, e.g.,
        "<gene1>,<gene2>,<gene3>". For weighted gene set, the 'GENESET' column contains
        gene names and weights, e.g., "<gene1>:<weight1>,<gene2>:<weight2>,<gene3>:<weight3>".
    src_species : str, default=None
        Source species, must be either 'mmusculus' or 'hsapiens' if not None
    dst_species : str, default=None
        Destination species, must be either 'mmusculus' or 'hsapiens' if not None
    to_intersect : List[str], default=None.
        Gene list to intersect with the input .gs file.

    Returns
    -------
    dict_gs : dict
        Dictionary of gene sets: {
            trait1: (gene_list, gene_weight_list),
            trait2: (gene_list, gene_weight_list),
            ...
        }
    """

    assert (src_species is None) == (
        dst_species is None
    ), "src_species and dst_species must be both None or not None"

    # Load homolog map dict_map; only needed when src_species and dst_species
    # are not None and different.
    if ((src_species is not None) & (dst_species is not None)) and (
        src_species != dst_species
    ):
        dict_map = load_homolog_mapping(src_species, dst_species)  # type: ignore
    else:
        dict_map = None  # type: ignore

    # Load gene set file
    dict_gs = {}
    df_gs = pd.read_csv(gs_path, sep="\t")
    for i, (trait, gs) in df_gs.iterrows():
        gs_info = [g.split(":") for g in gs.split(",")]
        if np.all([len(g) == 1 for g in gs_info]):
            # if all genes are weighted uniformly
            dict_weights = {g[0]: 1.0 for g in gs_info}
        elif np.all([len(g) == 2 for g in gs_info]):
            # if all genes are weighted by their weights
            dict_weights = {g[0]: float(g[1]) for g in gs_info}
        else:
            raise ValueError(f"gene set {trait} contains genes with invalid format")

        # Convert species if needed
        # convert gene list to homologs, if gene can not be mapped, remove it
        # in both gene list and gene weight
        if dict_map is not None:
            dict_weights = {
                dict_map[g]: w for g, w in dict_weights.items() if g in dict_map
            }

        # Intersect with other gene sets
        if to_intersect is not None:
            to_intersect = set(to_intersect)
            dict_weights = {g: w for g, w in dict_weights.items() if g in to_intersect}

        gene_list = list(dict_weights.keys())
        dict_gs[trait] = (
            gene_list,
            [dict_weights[g] for g in gene_list],
        )

    return dict_gs
