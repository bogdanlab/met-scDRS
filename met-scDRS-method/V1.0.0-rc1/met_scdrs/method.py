from scipy import sparse
import numpy as np
import time
from met_scdrs.util import get_memory

def preprocess(
    h5ad_obj,
    method : str = 'inverse',
    variance_clip : int = 5
):
    """
    Preprocess the methylation cell by gene data

    Parameters
    ----------
    h5ad_obj : anndata.AnnData
        AnnData object loaded from h5ad_file path
    method : str
        methodology to preprocess the single cell methylation matrix, supported methods:
        "inverse" : inverse the fraction into 1 - X
    variance_clip : int
        only genes with greater than specified percentile will be retained
        Default is 5 (i.e., remove low-variance genes under the 5th percentile).
    """
    # print header
    print(f'\n\n\n PREPROCESSING')
    initial_time = time.time()
    header = "--preprocess_method %s \\\n" % method
    header += "--variance_clip %s \\\n" % variance_clip
    print(header)
    
    # obtain the memory usage:
    print(f"Initiating preprocess, memory usage: {get_memory():.2f} MB")
    
    # define method (inverse): 
    def inverse_method(h5ad_obj):
        if sparse.issparse(h5ad_obj.X):
            print('h5ad object is recognized as sparse object')
            h5ad_obj.X.data = 1 - h5ad_obj.X.data
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
    methods = {'inverse' : inverse_method}
    preprocessed_data = methods[method](h5ad_obj)
    
    # method to usage 
    # print our usage information:
    print(f'Preprocess completed, elapsed time: {(time.time() - initial_time):.3f} seconds')
    print(f"Finished preprocess, memory usage: {get_memory():.2f} MB")
    return preprocessed_data


def score_cells(input_path):
    print(f"[met-scdrs] Input path received: {input_path}")