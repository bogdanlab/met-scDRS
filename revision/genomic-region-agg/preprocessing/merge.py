# conda activate allcools

import anndata as ad
import re
import os
import gc
import pandas as pd
from tqdm.auto import tqdm
import numpy as np
from datetime import date
today = date.today().strftime("%m%d%Y")

pd.set_option('display.max_columns', None)

# specify the tien files:
tien_10k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien10k/Tien10k_07092026.mcds"
tien_20k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien20k/Tien20k_07092026.mcds"
tien_30k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien30k/Tien30k_07092026.mcds"
tien_40k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien40k/Tien40k_07092026.mcds"
tien_50k = "/u/project/geschwind/lixinzhe/data/GSE215353/Tien50k/Tien50k_07092026.mcds"
files_to_loop = [tien_10k, tien_20k, tien_30k, tien_40k, tien_50k]
output_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted'
os.makedirs(f"{output_dir}/merged/", exist_ok = True)

# initiate the empty list:
file_path = []
file_base = []
file_mc_type = []
file_feature_space = []

# obtain which file information on what is the path, 10K, mc_type and feature space
for file in files_to_loop:
    for feature_space in ['promoter', 'exon', 'intron']:
        for mc_type in ['CHN', 'CGN']:
            output_base = re.sub('.mcds', '', re.sub('.*/', '', file))
            file_path.append(f"{output_dir}/{output_base}_{mc_type}_{feature_space}_extracted.h5ad")
            file_base.append(output_base)
            file_mc_type.append(mc_type)
            file_feature_space.append(feature_space)

# group this into a data frame:
file_info = pd.DataFrame({
    "file_path": file_path,
    "file_base": file_base,
    "file_mc_type": file_mc_type,
    "file_feature_space": file_feature_space
    }
)

###########################################################################################
######        Aggregate within each feature space, and mc type across base           ######
###########################################################################################
file_groups = file_info.groupby(["file_mc_type", "file_feature_space"])

# these variables are gene meta that are specific to each batch of the merge, so we should rename:
batch_var_columns = [
    "total_cov",
    "total_mc",
    "global_frac"
]

# for each group, aggregate:
for (mc_type, feature), sub in file_groups:
    
    # initiate empty adata:
    merged_adata = None
    
    # print a message:
    print(f'Merging within {mc_type} and {feature}')
    print('')
    
    for _, row in tqdm(sub.iterrows(), total = len(sub), desc = 'Merging cells'):
        current_path = row['file_path']
        current_base = row['file_base']
        
        current_adata = ad.read_h5ad(current_path)
        current_adata.obs["dataset"] = current_base
        current_adata.var = current_adata.var.rename(
            columns= {
                col: f'{col}_{current_base}'
                for col in batch_var_columns
                if col  in current_adata.var.columns    
            }
        )
        
        # check on the unique cell and gene name:
        assert current_adata.obs_names.is_unique, f"duplicated obs names observed in the {current_path}"
        assert current_adata.var_names.is_unique, f"duplicated var names observed in the {current_path}"
        
        if merged_adata is None:
            # first file is the merged data if none:
            merged_adata = current_adata
            del current_adata
            gc.collect()
        
        else:
            # perform a check on the variable name:
            same_features = set(merged_adata.var_names) == set(current_adata.var_names)
            same_order = merged_adata.var_names.equals(current_adata.var_names)
            if not same_features:
                print(f"Warning: feature sets differ in {current_base}")
            elif not same_order:
                print(f"Feature order differs in {current_base}, but names are the same")
                    
            # merge the accumulated data with the new file
            merged_adata = ad.concat(
                [merged_adata, current_adata],
                axis = 0,
                join = 'inner',
                merge="unique"
                )
            
            # perform sanity checks on unique names and same order features:
            assert merged_adata.obs_names.is_unique, 'cells name is duplicated between merges'
            
            # delete the unused objects:
            del current_adata
            gc.collect()
    
    # verify the current number of covariates are present per batch:
    cov_cols = [col for col in merged_adata.var.columns if col.startswith('total_cov_Tien')]
    mc_cols = [col for col in merged_adata.var.columns if col.startswith('total_mc_Tien')]
    
    assert len(cov_cols) == len(sub), (
        f"Expected {len(sub)} total_cov batch columns, "
        f"but found {len(cov_cols)}"
        )
    assert len(mc_cols) == len(sub), (
        f"Expected {len(sub)} total_mc batch columns, "
        f"but found {len(mc_cols)}"
        )
    
    # for the merged result, add in cov and mc:
    merged_adata.var['total_cov'] = (merged_adata.var[cov_cols].sum(axis = 1))
    merged_adata.var['total_mc'] = (merged_adata.var[mc_cols].sum(axis = 1))
    
    # calculate the merged results total percent methylation
    merged_adata.var["global_frac"] = np.divide(
        merged_adata.var["total_mc"].to_numpy(),
        merged_adata.var["total_cov"].to_numpy(),
        out=np.full(merged_adata.n_vars, np.nan),
        where=merged_adata.var["total_cov"].to_numpy() > 0
    )
    
    # check numerics:
    assert (merged_adata.var["total_cov"] >= 0).all()
    assert (merged_adata.var["total_mc"] >= 0).all()
    assert (merged_adata.var['total_mc'] <= merged_adata.var['total_cov']).all()
    
    # output:
    merged_adata.write_h5ad(f"{output_dir}/merged/merged_{today}_{mc_type}_{feature}_raw.h5ad")