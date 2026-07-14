# create meta data for the 50K data for visualization
import pandas as pd
import re

# load in the original meta data:
meta_path = '/u/project/cluo/heffel/BICAN3/REVISION/metadata_passQC_05212026.tsv.gz'
meta = pd.read_table(meta_path, sep = '\t', index_col = 0)

# load in one risk score:
risk_score_50k_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/baseline_promoter_CHN/PASS_MDD_Howard2019.score.gz'
risk_score_50k = pd.read_table(risk_score_50k_path, sep = '\t', index_col = 0)
risk_score_50k.index = [re.sub('.allc.tsv.gz', '', cell_name) for cell_name in risk_score_50k.index]

# check if the 50k index is in the meta:
risk_score_50k.index.isin(meta.index).sum() == len(risk_score_50k)
# 5 cells are missing from the meta

# output the subset of meta data:
common_index = list(set(risk_score_50k.index).intersection(set(meta.index)))
meta_50k = meta.loc[common_index, :]
meta_50k.to_csv('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv', sep = '\t')

