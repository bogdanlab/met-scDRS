import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date

# specify path to read:
input_dir = '/u/project/geschwind/lixinzhe/data/GSE215353/extracted/merged_QCed'
today = date.today().strftime("%Y-%m-%d")
plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'


meta = pd.read_table('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv', sep = '\t', index_col = 0)

def plot_mc_cov_density(
    adata_meta,
    mc_column: str = "total_mc",
    cov_column: str = "total_cov",
):
    fig, ax = plt.subplots(figsize=(8, 5))

    sns.kdeplot(
        data=adata_meta,
        x=mc_column,
        label="Methylated count",
        ax=ax,
    )

    sns.kdeplot(
        data=adata_meta,
        x=cov_column,
        label="Coverage count",
        ax=ax,
    )

    ax.set_xlabel("Count")
    ax.set_ylabel("Density")
    ax.legend()

    fig.tight_layout()

    return fig, ax

def plot_methyl_density(
    adata_meta,
    cell_type_column = 'cell_type',
    mc_column: str = "total_mc",
    cov_column: str = "total_cov"
    ):
    plot_df = adata_meta.copy()
    
    # compute the global fraction
    plot_df["methylation_fraction"] = (
        plot_df[mc_column] / plot_df[cov_column]
    )
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # draw a fraction density plot
    sns.kdeplot(
        data=plot_df,
        x="methylation_fraction",
        hue=cell_type_column,
        common_norm=False,
        ax=ax,
    )

    ax.set_xlabel("Methylation fraction")
    ax.set_ylabel("Density")
    fig.tight_layout()
    return fig, ax


features = ['promoter', 'exon', 'intron']
mc_types = ['CHN', 'CGN']

for mc_type in mc_types:
    for feature in features:
        h5ad_path = f"{input_dir}/merged_071526_{mc_type}_{feature}_QC.h5ad"
        current_h5ad = ad.read_h5ad(h5ad_path)
        current_h5ad.obs.index = [re.sub(r"\.allc\.tsv\.gz$", '', cell_name) for cell_name in current_h5ad.obs.index]
        common_index = list(set(current_h5ad.obs.index).intersection(set(meta.index)))

        current_adata_obs = current_h5ad.obs.copy()
        current_adata_obs = current_adata_obs.loc[common_index, :]
        current_adata_obs['cell_type'] = meta.loc[common_index, 'L1']
        
        # make a plot of methylated fraction separated by cell class:
        fig, ax = plot_methyl_density(current_adata_obs)
        ax.set_title(f"{feature}\n{mc_type}")

        fig.savefig(
            f"{plot_dir}/{today}-methylation_density_by_cell_class-{mc_type}-{feature}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)
        
        # make another plot of methylated coverage and counts:
        fig, ax = plot_mc_cov_density(current_adata_obs)
        ax.set_title(f"{feature}\n{mc_type}")

        fig.savefig(
            f"{plot_dir}/{today}-methylation_counts-{mc_type}-{feature}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)