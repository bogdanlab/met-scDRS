###########################################################################################
######                                    PREAMBLE                                   ######
###########################################################################################
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from datetime import date
today = date.today()
from statsmodels.stats.multitest import multipletests

# load in the scores:
# we can use MDD as an example
disease_score_name = 'PASS_MDD_Howard2019.score.gz'

# path for gene body, after QC but no correction
gene_body_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/'

# also for different aggregation after QC but no correction:
promoter_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CHN/'
intron_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CHN/'
exon_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CHN/'

# load in the score
gene_body_mdd = pd.read_table(f"{gene_body_path}{disease_score_name}", index_col = 0)
promoter_mdd = pd.read_table(f"{promoter_path}{disease_score_name}", index_col = 0)
exon_mdd = pd.read_table(f"{exon_path}{disease_score_name}", index_col = 0)
intron_mdd = pd.read_table(f"{intron_path}{disease_score_name}", index_col = 0)

###########################################################################################
######                                    Visualization                              ######
###########################################################################################
def visualize_XY(
    X,
    Y,
    x_name = "Gene body",
    y_name = 'Y',
    score_col = 'zscore',
    ax = None,
    alpha = 0.35,
    s = 8,
    color = None
    ):
    """
    X-Y scatter plot comparing between other region vs gene body
    """
    # subst to the common cells and the score column:
    common_cells = X.index.intersection(Y.index)
    plot_df = pd.DataFrame({
        "x": X.loc[common_cells, score_col],
        "y": Y.loc[common_cells, score_col],
    }).dropna()
    
    # inherit plotting axis if needed:
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    else:
        fig = ax.figure
    
    # get the pearson and spearman correlation:
    r, p = pearsonr(plot_df["x"], plot_df["y"])
    rho, sp = spearmanr(plot_df["x"], plot_df["y"])
    
    # make the scatter plot:
    ax.scatter(plot_df["x"], plot_df["y"], alpha=alpha, s=s, color = color)
    ax.axhline(0, color="gray", lw=0.8, alpha=0.5)
    ax.axvline(0, color="gray", lw=0.8, alpha=0.5)
    
    ax.set_xlabel(f"{x_name} zscore")
    ax.set_ylabel(f"{y_name} zscore")
    ax.set_title(
        f"{x_name} vs {y_name}\n"
        f"n={len(plot_df):,}, r={r:.3f}, rho={rho:.3f}"
    )
    return fig, ax


# use the defined function to loop over the y and x
y_dict = {
    "Promoter": promoter_mdd,
    "Exon": exon_mdd,
    "Intron": intron_mdd,
}
fig, axes = plt.subplots(1, len(y_dict), figsize=(15, 4), constrained_layout=True)

for ax, (y_name, Y) in zip(axes, y_dict.items()):
    visualize_XY(
        X=gene_body_mdd,
        Y=Y,
        x_name="Gene body",
        y_name=y_name,
        ax=ax,
        color="#F58518",
    )

plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'
fig.savefig(
    f"{plot_dir}{today}_nonCpG_MDD_zscore_gene_body_vs_promoter_exon_intron.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

###########################################################################################
######                             Also make a bar plot                              ######
###########################################################################################
def count_sig_cells(
    df,
    p_col="pval",
    alpha=0.05,
    method="fdr_bh",
):
    p_adj = multipletests(df[p_col].dropna(), method=method)[1]
    return (p_adj < alpha).sum()

score_dict = {
    "Gene body": gene_body_mdd,
    "Promoter": promoter_mdd,
    "Exon": exon_mdd,
    "Intron": intron_mdd,
}

sig_counts = pd.Series({
    name: count_sig_cells(df, p_col="pval", alpha=0.1, method="fdr_bh")
    for name, df in score_dict.items()
})
sig_counts = sig_counts.sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(5, 4))

ax.bar(sig_counts.index, sig_counts.values, color="#F58518")
ax.set_ylabel("Number of significant cells")
ax.set_xlabel("Region")
ax.set_title("Significant MDD-scoring cells after FDR correction")

for i, value in enumerate(sig_counts.values):
    ax.text(i, value, f"{value:,}", ha="center", va="bottom")

fig.savefig(
    f"{plot_dir}{today}_nonCpG_MDD_significant_cells_after_FDR.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

###########################################################################################
######                                    Extend to CpG                              ######
###########################################################################################
# path for gene body, after QC but no correction
gene_body_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/'

# also for different aggregation after QC but no correction:
promoter_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_promoter_CGN/'
intron_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_intron_CGN/'
exon_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_exon_CGN/'

# load in the score
gene_body_mdd = pd.read_table(f"{gene_body_path}{disease_score_name}", index_col = 0)
promoter_mdd = pd.read_table(f"{promoter_path}{disease_score_name}", index_col = 0)
exon_mdd = pd.read_table(f"{exon_path}{disease_score_name}", index_col = 0)
intron_mdd = pd.read_table(f"{intron_path}{disease_score_name}", index_col = 0)

# use the defined function to loop over the y and x
y_dict = {
    "Promoter": promoter_mdd,
    "Exon": exon_mdd,
    "Intron": intron_mdd,
}
fig, axes = plt.subplots(1, len(y_dict), figsize=(15, 4), constrained_layout=True)

for ax, (y_name, Y) in zip(axes, y_dict.items()):
    visualize_XY(
        X=gene_body_mdd,
        Y=Y,
        x_name="Gene body",
        y_name=y_name,
        ax=ax,
        color="#4C78A8",
    )

plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'
fig.savefig(
    f"{plot_dir}{today}_MDD_CpG_zscore_gene_body_vs_promoter_exon_intron.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

score_dict = {
    "Gene body": gene_body_mdd,
    "Promoter": promoter_mdd,
    "Exon": exon_mdd,
    "Intron": intron_mdd,
}

sig_counts = pd.Series({
    name: count_sig_cells(df, p_col="pval", alpha=0.1, method="fdr_bh")
    for name, df in score_dict.items()
})
sig_counts = sig_counts.sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(5, 4))

ax.bar(sig_counts.index, sig_counts.values, color = "#4C78A8")
ax.set_ylabel("Number of significant cells")
ax.set_xlabel("Region")
ax.set_title("Significant MDD-scoring cells after FDR correction")

for i, value in enumerate(sig_counts.values):
    ax.text(i, value, f"{value:,}", ha="center", va="bottom")

fig.savefig(
    f"{plot_dir}{today}_MDD_CpG_significant_cells_after_FDR.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)
