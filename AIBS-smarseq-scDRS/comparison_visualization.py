# the goal of the script is to compare scDRS with met-scDRS
# load in some packages
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
import anndata as ad
import numpy as np

pd.set_option('display.max_columns', None)
year_month_day = date.today().strftime("%Y-%m-%d")

###########################################################################################
######                                    DATA LOADING                               ######
###########################################################################################
# load in the results for the scDRS:
trx_dir = "/u/project/geschwind/lixinzhe/scDRS-output/scDRS-output/AIBS-psych-trait-scDRS/without_cov/"
# trx_dir = '/u/project/geschwind/lixinzhe/scDRS-output/scDRS-output/AIBS-psych-trait-scDRS/with_cov/'
met_dir = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_gene_body_CHN/'
plot_dir = '/u/home/l/lixinzhe/project-geschwind/plot/'

# load in five different diseases:
traits = [
    'PASS_Alzheimers_Jansen2019.score.gz',
    'PASS_BIP_Mullins2021.score.gz',
    'PASS_MDD_Howard2019.score.gz',
    'PASS_Schizophrenia_Pardinas2018.score.gz',
    'UKB_460K.body_HEIGHTz.score.gz'
]

trx_score = {}
met_score = {}

for trait in traits:
    # specify path
    trx_path = f"{trx_dir}{trait}"
    met_path = f"{met_dir}{trait}"
    
    # load in:
    trait_name = trait.replace('.score.gz', '')
    trx_score[trait_name] = pd.read_table(trx_path, sep = '\t', index_col = 0)
    met_score[trait_name] = pd.read_table(met_path, sep = '\t', index_col = 0)
    
    # Benjamini-Hochberg FDR correction within this trait and modality
    trx_score[trait_name]["fdr"] = multipletests(trx_score[trait_name]["pval"], method="fdr_bh")[1]
    met_score[trait_name]['fdr'] = multipletests(met_score[trait_name]["pval"], method="fdr_bh")[1]
    
    # change the index name:
    met_score[trait_name].index = met_score[trait_name].index.str.replace('.allc.tsv.gz', '')

# also load in the meta data for both dataset:
meta_trx = pd.read_table('/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/metadata.csv', sep = ',', index_col = 0)
meta_methyl = pd.read_table('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv', sep = '\t', index_col = 0)
fdr_cutoff = 0.1


###########################################################################################
######                                    Comparison                                 ######
###########################################################################################
# convert the cell class label so that they are more comparable:
trx_class_map = {
    "GABAergic": "Inhibitory",
    "Glutamatergic": "Excitatory",
    "Non-neuronal": "Non-neuronal",
}

methyl_class_map = {
    "Inh": "Inhibitory",
    "Exc": "Excitatory",
    "Glial": "Non-neuronal",
    "NN": "Non-neuronal",
}

meta_trx["broad_cell_class"] = meta_trx["class_label"].map(trx_class_map)
meta_methyl["broad_cell_class"] = meta_methyl["L1"].map(methyl_class_map)
print(meta_trx["broad_cell_class"].value_counts(dropna=False))
print(meta_methyl["broad_cell_class"].value_counts(dropna=False))

### Compare the proportion of cell class that are significant in 5 different traits
meta_info = {
    "scDRS": {
        "score": trx_score,
        "meta": meta_trx,
    },
    "met-scDRS": {
        "score": met_score,
        "meta": meta_methyl,
    },
}

prop_sig_rows = []
for trait_name in trx_score:
    for modality, info in meta_info.items():
        score_df = info["score"][trait_name].copy()
        meta_df = info["meta"]

        df = score_df[["pval", "fdr"]].join(
            meta_df[["broad_cell_class"]],
            how="inner"
        )
        
        df = df.dropna(subset=["broad_cell_class", "pval"])
        
        # get the significant label
        df["significant"] = df["fdr"] < fdr_cutoff
        
        # for each of the broad cell class, summaize by number of significant cells:
        summary = (
            df.groupby("broad_cell_class")
            .agg(
                n_cells=("pval", "size"),
                n_significant=("significant", "sum"),
                prop_significant=("significant", "mean"),
                percent_significant=("significant", "mean"),
                min_pval=("pval", "min"),
                min_fdr=("fdr", "min"),
            )
            .reset_index()
        )
        summary["percent_significant"] = summary["percent_significant"] * 100
        summary["trait"] = trait_name
        summary["modality"] = modality
        summary["fdr_cutoff"] = fdr_cutoff
        
        # append the summary:
        prop_sig_rows.append(summary)

# get the proportion:
prop_sig_by_class = pd.concat(prop_sig_rows, ignore_index=True)

###########################################################################################
######                                    Visualization                              ######
###########################################################################################
class_order = ["Excitatory", "Inhibitory", "Non-neuronal"]

g = sns.catplot(
    data=prop_sig_by_class,
    kind="bar",
    x="broad_cell_class",
    y="percent_significant",
    hue="modality",
    col="trait",
    col_wrap=3,
    order=class_order,
    height=4,
    aspect=1.2,
    errorbar=None,
)

g.set_axis_labels("", "% significant cells, FDR < 0.1")
g.set_titles("{col_name}")

for ax in g.axes.flat:
    ax.tick_params(axis="x", rotation=30)

plt.subplots_adjust(top=0.88)
g.fig.suptitle("Proportion of FDR-significant cells by broad cell class")
out_path_png = os.path.join(plot_dir, f"{year_month_day}_prop_fdr_significant_cells_by_class.png")
g.fig.savefig(out_path_png, dpi=300, bbox_inches="tight")
plt.close(g.fig)

###########################################################################################
######                               Visualize within non-neuronal                   ######
###########################################################################################
# load in the anndata
tsne_trx = pd.read_table("/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/tsne.csv", sep = ',', index_col = 0)
umap_met = pd.read_table("/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv", sep = ',', index_col = 0)

# specify traits:
traits_to_plot = {
    "Alzheimers": "PASS_Alzheimers_Jansen2019",
    "MDD": "PASS_MDD_Howard2019",
}

def plot_score_embedding(
    coord_df,
    score_df,
    x_col, 
    y_col,
    score_col, 
    title,
    out_prefix,
    vmin = None,
    vmax = None,
    cmap = 'RdBu_r',
):
    # subset to the common set of cells:
    shared_cells = coord_df.index.intersection(score_df.index)
    
    # specify the plot_df:
    plot_df = coord_df.loc[shared_cells, [x_col, y_col]].join(
        score_df.loc[shared_cells, [score_col]],
        how="inner",
    )
    plot_df = plot_df.dropna(subset=[x_col, y_col, score_col])
    
    # short circuit if there is no shared cells:
    if plot_df.shape[0] == 0:
        print(f"Skipping {title}: no shared cells")
        return
    
    # If the we never provided the vmin and vmax:
    if vmin is None or vmax is None:
        vmax = np.nanpercentile(np.abs(plot_df[score_col]), 99)
        vmin = -vmax
    
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        plot_df[x_col],
        plot_df[y_col],
        c=plot_df[score_col],
        s=3,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )
    
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_xticks([])
    ax.set_yticks([])
    
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(score_col)

    fig.tight_layout()

    png_path = os.path.join(plot_dir, f"{year_month_day}-{out_prefix}.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"{title}: n cells plotted = {plot_df.shape[0]}")
    print(f"Saved PNG: {png_path}")


# first get the global color bar:
score_col = "zscore"
all_scores = []

# concatenate across traits
for trait_label, trait_name in traits_to_plot.items():
    trx_shared = tsne_trx.index.intersection(trx_score[trait_name].index)
    met_shared = umap_met.index.intersection(met_score[trait_name].index)

    all_scores.append(trx_score[trait_name].loc[trx_shared, score_col])
    all_scores.append(met_score[trait_name].loc[met_shared, score_col])
all_scores = pd.concat(all_scores).dropna()

# get teh global high and low percentile
global_vmax = np.nanpercentile(np.abs(all_scores), 99)
global_vmin = -global_vmax

print("Global color scale:")
print("vmin:", global_vmin)
print("vmax:", global_vmax)


for trait_label, trait_name in traits_to_plot.items():
    # scDRS on transcriptomic tSNE
    plot_score_embedding(
        coord_df=tsne_trx,
        score_df=trx_score[trait_name],
        x_col="tsne_1",
        y_col="tsne_2",
        score_col="zscore",
        vmin = global_vmin,
        vmax = global_vmax,
        title=f"scDRS {trait_label} z-score on transcriptomic data",
        out_prefix=f"scdrs_{trait_label}_zscore_on_trx_tsne",
    )

    # met-scDRS on methylation UMAP
    plot_score_embedding(
        coord_df=umap_met,
        score_df=met_score[trait_name],
        x_col="X_umap",
        y_col="Y_umap",
        vmin = global_vmin,
        vmax = global_vmax,
        score_col="zscore",
        title=f"met-scDRS {trait_label} z-score on methylation data",
        out_prefix=f"met_scdrs_{trait_label}_zscore_on_met_umap",
    )

###########################################################################################
######              CELL CLASS EMBEDDINGS WITH SELECTED SUBTYPE ANNOTATIONS          ######
###########################################################################################

import matplotlib.patheffects as pe
cell_class_palette = {
    "Excitatory": "#4C78A8",
    "Inhibitory": "#F58518",
    "Non-neuronal": "#54A24B",
}

def plot_label_embedding_with_subtype_text(
    coord_df,
    meta_df,
    x_col,
    y_col,
    label_col,
    subtype_col,
    annotate_class,
    title,
    out_prefix,
    palette=cell_class_palette,
    min_cells_for_label=20,
):
    shared_cells = coord_df.index.intersection(meta_df.index)

    plot_df = coord_df.loc[shared_cells, [x_col, y_col]].join(
        meta_df.loc[shared_cells, [label_col, subtype_col]],
        how="inner",
    )

    plot_df = plot_df.dropna(subset=[x_col, y_col, label_col, subtype_col])

    if plot_df.shape[0] == 0:
        print(f"Skipping {title}: no shared cells")
        return

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    # Keep broad class colors unchanged
    for label, sub_df in plot_df.groupby(label_col):
        ax.scatter(
            sub_df[x_col],
            sub_df[y_col],
            s=3,
            color=palette.get(label, "lightgray"),
            label=label,
            linewidths=0,
            rasterized=True,
        )

    # Add subtype text only within the requested broad class
    anno_df = plot_df[plot_df[label_col] == annotate_class].copy()

    subtype_centers = (
        anno_df.groupby(subtype_col)
        .agg(
            x=(x_col, "median"),
            y=(y_col, "median"),
            n=(subtype_col, "size"),
        )
        .reset_index()
    )

    subtype_centers = subtype_centers[
        subtype_centers["n"] >= min_cells_for_label
    ]

    for _, row in subtype_centers.iterrows():
        ax.text(
            row["x"],
            row["y"],
            str(row[subtype_col]),
            fontsize=7,
            ha="center",
            va="center",
            color="black",
            path_effects=[
                pe.withStroke(linewidth=2.5, foreground="white")
            ],
        )

    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.legend(
        title=label_col,
        markerscale=4,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    fig.tight_layout()

    png_path = os.path.join(plot_dir, f"{year_month_day}-{out_prefix}.png")
    pdf_path = os.path.join(plot_dir, f"{year_month_day}-{out_prefix}.pdf")

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"{title}: n cells plotted = {plot_df.shape[0]}")
    print(f"Annotated {subtype_centers.shape[0]} subtype labels")
    print(f"Saved PNG: {png_path}")
    print(f"Saved PDF: {pdf_path}")

plot_label_embedding_with_subtype_text(
    coord_df=tsne_trx,
    meta_df=meta_trx,
    x_col="tsne_1",
    y_col="tsne_2",
    label_col="broad_cell_class",
    subtype_col="subclass_label",
    annotate_class="Non-neuronal",
    title="Transcriptomic tSNE colored by broad cell class; non-neuronal subtypes labeled",
    out_prefix="trx_tsne_broad_cell_class_non_neuronal_subclass_labels",
)

plot_label_embedding_with_subtype_text(
    coord_df=umap_met,
    meta_df=meta_methyl,
    x_col="X_umap",
    y_col="Y_umap",
    label_col="broad_cell_class",
    subtype_col="L3",
    annotate_class="Excitatory",
    title="Methylation UMAP colored by broad cell class; excitatory subtypes labeled",
    out_prefix="met_umap_broad_cell_class_excitatory_L3_labels",
)