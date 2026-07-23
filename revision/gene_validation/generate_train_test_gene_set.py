###########################################################################################
######                               Generate train test gene set                    ######
###########################################################################################
# load in packages:
import anndata as ad
import pickle
import pandas as pd
import os
import numpy as np
import gc

pd.set_option('display.max_columns', None)

# load in adata:
with open("/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/intermediate_merged_071726_CHN_gene_body_QC.pkl", "rb") as f:
    data = pickle.load(f)

out_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/preprocessed_pkl_gene_batch/'
os.makedirs(out_dir, exist_ok = True)
os.makedirs(f"{out_dir}/gene_batch/", exist_ok = True)

###########################################################################################
######                                    get the genes                              ######
###########################################################################################
# get the genes from the original data:
genes = np.array(data.var_names)

# get the rng:
rng = np.random.default_rng(0)
shuffled_genes = rng.permutation(genes)

# split the shuffled genes into 10 different folds
heldout_folds = np.array_split(shuffled_genes, 10)

for fold_idx, heldout_genes in enumerate(heldout_folds, start=1):
    print(f"Processing fold {fold_idx}")
    keep_genes = genes[~np.isin(genes, heldout_genes)]
    
    # save gene list:
    pd.Series(heldout_genes).to_csv(f'{out_dir}/gene_batch/fold{fold_idx}_held_out_gene_list.csv', index = False, header = False)
    pd.Series(keep_genes).to_csv(f'{out_dir}/gene_batch/fold{fold_idx}_keep_gene_list.csv', index = False, header = False)
    
    # make the 90% gene dataset:
    subset = data[:, keep_genes].copy()
    
    # make sure scDRS gene stats only include kept genes:
    gene_stats = subset.uns["SCDRS_PARAM"]["GENE_STATS"]
    subset.uns["SCDRS_PARAM"]["GENE_STATS"] = gene_stats.loc[keep_genes]
    
    # add checks:
    assert subset.n_vars == len(keep_genes)
    assert subset.uns["SCDRS_PARAM"]["GENE_STATS"].shape[0] == len(keep_genes)
    assert np.array_equal(subset.var_names.to_numpy(), keep_genes)
    assert np.array_equal(
        subset.uns["SCDRS_PARAM"]["GENE_STATS"].index.to_numpy(),
        keep_genes
    )
    
    # Save the actual AnnData object as pickle
    with open(f"{out_dir}/fold_{fold_idx}_data.pkl", "wb") as f:
        pickle.dump(subset, f)
    
    # free up memory:
    del subset
    gc.collect()

# concatenate the held out to add checkes:
all_heldout = np.concatenate(heldout_folds)

assert len(all_heldout) == len(genes)
assert len(np.unique(all_heldout)) == len(genes)
assert set(all_heldout) == set(genes)
print("All checks passed.")
    