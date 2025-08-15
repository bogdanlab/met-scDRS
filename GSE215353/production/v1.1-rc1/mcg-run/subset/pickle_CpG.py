import scanpy as sc
import numpy as np
import pickle
import met_scdrs

h5ad = "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/2025-08-05-mcg-subset-knn-aggregate.h5ad"
data = sc.read_h5ad(h5ad)

# process the result:
met_scdrs.preprocess(
    data,
    cov=None,
    n_mean_bin=10,
    n_var_bin=10,
    n_length_bin = 10,
    copy=False,
    weight_option="inv_std",
    ctrl_match_key="mean_var_length",
    verbose = True)

# pickle the result:
output_pickle = "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/2025-08-05-mcg-subset-knn-aggregate.pkl"
with open(output_pickle, "wb") as f: pickle.dump(data, f)