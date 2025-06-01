from pathlib import Path
import scanpy as sc
import anndata
import numpy as np

def convert_species_name(species):
    if species in ["Mouse", "mouse", "Mus_musculus", "mus_musculus", "mmusculus"]:
        return "mmusculus"
    if species in ["Human", "human", "Homo_sapiens", "homo_sapiens", "hsapiens"]:
        return "hsapiens"
    raise ValueError("species name '%s' is not supported" % species)

def check_folder(folder_path):
    folder = Path(folder_path)
    if folder.exists() and folder.is_dir():
        return True
    else:
        return False

def load_h5ad(h5ad_file : str):
    """
    Load in h5ad file and check if there is NA in the data:

    Parameters
    ----------
    h5ad_file : str
        Path to the single-cell `.h5ad` file. The `.X` attribute should contain a cell-by-gene 
        methylation matrix.
    """
    adata = anndata.read_h5ad(h5ad_file)
    
    # check inputs (1) no NaN in adata.X
    if np.isnan(adata.X.sum()):
        raise ValueError(
            "h5ad expression matrix should not contain NaN. Please impute them beforehand."
        )
    return adata
