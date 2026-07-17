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

