# Purpose: compare a with b
# where a = QCed, mean var density, full regression on coverage metrics
# and b = QCed, mean var length, global methylation fraction regression

# make scatter plot to look at the risk score correlations

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

# path for gene body, after QC but no correction, mean var length
mean_var_length_no_regress_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/'

# path for gene body after QC, with rowSum correction, mean var length:
mean_var_length_rowsum_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/cov/'

# path for gene body, after QC, with full correction, mean var density:
mean_var_density_coverage_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CHN/full/'

# load in the score
mean_var_length_no_regress = pd.read_table(f"{mean_var_length_no_regress_path}{disease_score_name}", index_col = 0)
mean_var_length_rowsum = pd.read_table(f"{mean_var_length_rowsum_path}{disease_score_name}", index_col = 0)
mean_var_density_coverage = pd.read_table(f"{mean_var_density_coverage_path}{disease_score_name}", index_col = 0)

# define a function:
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
    
    # y = x reference line
    lim_min = min(plot_df["x"].min(), plot_df["y"].min())
    lim_max = max(plot_df["x"].max(), plot_df["y"].max())
    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        color="black",
        lw=1,
        linestyle="--",
        alpha=0.7,
        label="y = x",
    )
    
    ax.set_xlabel(f"{x_name} zscore")
    ax.set_ylabel(f"{y_name} zscore")
    ax.set_title(
        f"{x_name} vs {y_name}\n"
        f"n={len(plot_df):,}, r={r:.3f}, rho={rho:.3f}"
    )
    return fig, ax

# use the defined function to loop over the y and x
y_dict = {
    "Coverage_aware_correction": mean_var_density_coverage,
    "Coverage_unaware_correction": mean_var_length_rowsum
}
fig, axes = plt.subplots(1, len(y_dict), figsize=(15, 4), constrained_layout=True)

for ax, (y_name, Y) in zip(axes, y_dict.items()):
    visualize_XY(
        X=mean_var_length_no_regress,
        Y=Y,
        x_name="Baseline",
        y_name=y_name,
        ax=ax,
        color="#F58518",
    )

plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'
fig.savefig(
    f"{plot_dir}{today}_nonCpG_MDD_zscore_preprocessing_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

###########################################################################################
######                        Also Compare Non-CpG methylation                       ######
###########################################################################################
# load in the scores:
# we can use MDD as an example
disease_score_name = 'PASS_MDD_Howard2019.score.gz'

# path for gene body, after QC but no correction, mean var length
mean_var_length_no_regress_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/'

# path for gene body after QC, with rowSum correction, mean var length:
mean_var_length_rowsum_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CGN/cov/'

# path for gene body, after QC, with full correction, mean var density:
mean_var_density_coverage_path = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CGN/full/'

# load in the score
mean_var_length_no_regress = pd.read_table(f"{mean_var_length_no_regress_path}{disease_score_name}", index_col = 0)
mean_var_length_rowsum = pd.read_table(f"{mean_var_length_rowsum_path}{disease_score_name}", index_col = 0)
mean_var_density_coverage = pd.read_table(f"{mean_var_density_coverage_path}{disease_score_name}", index_col = 0)

# use the defined function to loop over the y and x
y_dict = {
    "Coverage_aware_correction": mean_var_density_coverage,
    "Coverage_unaware_correction": mean_var_length_rowsum
}
fig, axes = plt.subplots(1, len(y_dict), figsize=(15, 4), constrained_layout=True)

for ax, (y_name, Y) in zip(axes, y_dict.items()):
    visualize_XY(
        X=mean_var_length_no_regress,
        Y=Y,
        x_name="Baseline",
        y_name=y_name,
        ax=ax,
        color="#4C78A8",
    )

plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'
fig.savefig(
    f"{plot_dir}{today}_CpG_MDD_zscore_preprocessing_comparison.png",
    dpi=300,
    bbox_inches="tight"
)