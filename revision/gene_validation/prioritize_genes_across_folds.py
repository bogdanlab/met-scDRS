###########################################################################################
######                                 PREAMBLE                                      ######
###########################################################################################
# load in packages:
import pickle
import pandas as pd
import numpy as np
import anndata as ad
import re
import os
pd.set_option('display.max_columns', None)
from scipy.stats import pearsonr
from tqdm.auto import tqdm

# load in adata:
with open("/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/intermediate_merged_071726_CHN_gene_body_QC.pkl", "rb") as f:
    full_methylation = pickle.load(f)
meta = pd.read_table('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv', sep = '\t', index_col = 0)

###########################################################################################
######                                    Define function                            ######
###########################################################################################
def get_correlation(adata, score, test_genes, gs_file, trait):
    # assert the score index is in the adata:
    assert (score.index.isin(adata.obs_names)).all()
    
    # intersect betweeen the test genes and the gs file within the trait:
    gene_set_in_trait = (
        gs_file.loc[gs_file.TRAIT == trait, "GENESET"]
        .dropna()
        .str.split(",")
        .explode()
        .str.split(":")
        .str[0]
        .tolist()
    )
    
    # get the genes that we actually need correlation to:
    test_genes_set = set(test_genes)
    genes_in_scope = [g for g in gene_set_in_trait if g in test_genes_set and g in adata.var_names]
    
    if len(genes_in_scope) == 0:
        return pd.DataFrame(columns=["gene", "corr"])
    # subset the data:
    working_adata = adata[score.index, genes_in_scope]
    X = working_adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    y = score.loc[working_adata.obs_names, "zscore"].to_numpy(dtype = float)
    
    # compute correlation:
    corr_vec = np.corrcoef(X, y, rowvar=False)[-1, :-1]
    # # check if they are close:
    # # for loop version
    # corr_loop = np.array([
    #     pearsonr(X[:, i], y).statistic
    #     for i in range(X.shape[1])
    # ])
    # np.allclose(corr_vec, corr_loop, equal_nan=True)
    
    # return the result:
    corr_result = pd.DataFrame({
        "gene": working_adata.var_names,
        "corr": corr_vec,
    })
    return corr_result


###########################################################################################
######                                    Load in data                               ######
###########################################################################################
batch_names = [f'fold_{fold_index}' for fold_index in range(1, 11)]
gene_train_test_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/preprocessed_pkl_gene_batch/gene_batch/'

# load in the gs files:
gs_path = "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs"
gs_file = pd.read_table(gs_path, sep = '\t')

# for each of the batch, get the kept genes and the fold index:
all_cor = []
for batch_name in tqdm(batch_names, desc = 'Processing across folds'):
    # load in the genes in test batch:
    fold_index = re.sub('_', '', batch_name)
    train_gene_path = f"{gene_train_test_dir}{fold_index}_keep_gene_list.csv"
    test_gene_path = f"{gene_train_test_dir}{fold_index}_held_out_gene_list.csv"

    # load in the genes:
    train_gene = pd.read_table(train_gene_path, header = None).iloc[:, 0].tolist()
    test_gene = pd.read_table(test_gene_path, header = None).iloc[:, 0].tolist()
    
    # set the directory name:
    load_dir = f'/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/batch/{batch_name}/'

    # load in the score:
    existing_scores = [score_file for score_file in os.listdir(load_dir) if score_file.endswith(".score.gz")]
    for score in tqdm(existing_scores, desc = 'Processing across scores'):
        trait = score.replace(".score.gz", "")
        score_path = f'{load_dir}{score}'
        load_score = pd.read_table(score_path, sep = '\t', index_col = 0)
    
        # compute the correlation between the score and the methylation on the test gene:
        cor_in_fold = get_correlation(adata = full_methylation, score = load_score, test_genes = test_gene, gs_file = gs_file, trait = trait)
        cor_in_fold["batch"] = batch_name
        cor_in_fold["trait"] = trait
        all_cor.append(cor_in_fold)

# append all the correlations:
if len(all_cor) > 0:
    all_cor = pd.concat(all_cor, ignore_index=True)
else:
    all_cor = pd.DataFrame(columns=["gene", "corr", "batch", "trait"])

all_cor.to_csv(
    "/u/home/l/lixinzhe/project-geschwind/port/scDRS/test_gene_methylation_score_correlations.tsv",
    sep="\t",
    index=False,
)
