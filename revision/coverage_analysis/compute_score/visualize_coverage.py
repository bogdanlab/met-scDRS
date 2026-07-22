import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
today = date.today().strftime("%Y-%m-%d")
plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'

mcg_cov = pd.read_table(
    "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CGN_gene_body_full.cov",
    sep = '\t',
    index_col = 0
    )

mch_cov = pd.read_table(
    "/u/project/geschwind/lixinzhe/data/GSE215353/extracted_gene_body/merged_QCed/CHN_gene_body_full.cov",
    sep = '\t',
    index_col = 0
    )

meta_50k = pd.read_table(
    '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv',
    sep = '\t',
    index_col = 0
    )

mcg_cov.index = [re.sub('.allc.tsv.gz', '', index_id) for index_id in mcg_cov.index]
mch_cov.index = [re.sub('.allc.tsv.gz', '', index_id) for index_id in mch_cov.index]

# match them by ID:
common_cells = (
    mcg_cov.index
    .intersection(mch_cov.index)
    .intersection(meta_50k.index)
)

# Preserve metadata order
common_cells = meta_50k.index[meta_50k.index.isin(common_cells)]

mcg_cov = mcg_cov.loc[common_cells].copy()
mch_cov = mch_cov.loc[common_cells].copy()
meta_50k = meta_50k.loc[common_cells].copy()

assert mch_cov.index.equals(mcg_cov.index)
assert mch_cov.index.equals(meta_50k.index)

# add cell class
mcg_cov["L1"] = meta_50k["L1"]
mch_cov["L1"] = meta_50k["L1"]

###########################################################################################
######                                    Plot                                       ######
###########################################################################################
plt.figure(figsize=(8, 5))

sns.kdeplot(
    data=mcg_cov,
    x="coverage",
    hue="L1",
    common_norm=False,
)

plt.xlabel("Total CG coverage")
plt.title("CG coverage by cell class")
plt.tight_layout()
plt.savefig(
    f"{plot_dir}/{today}_mcg_gene_body_coverage_by_L1.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
plt.close()



plt.figure(figsize=(8, 5))

sns.kdeplot(
    data=mch_cov,
    x="coverage",
    hue="L1",
    common_norm=False,
)

plt.xlabel("Total CH coverage")
plt.title("CH coverage by cell class")
plt.tight_layout()
plt.savefig(
    f"{plot_dir}/{today}_mch_gene_body_coverage_by_L1.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()