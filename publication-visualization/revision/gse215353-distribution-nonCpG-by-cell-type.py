# purpose: draw supplementary figures, get the distribution between non-CpG post process with respect to cell type

# load in data:
# load in packages:
import pickle
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
today = datetime.today().strftime("%Y-%m-%d")

###########################################################################################
######                                   LOADING DATA                                ######
###########################################################################################
# load in the meta data:
meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')

# load in the actual data:
intermediate_pkl = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.pkl'
with open(intermediate_pkl, 'rb') as f:
    processed_h5ad = pickle.load(f)

# use polars to output the file:
h5ad_pl = pl.DataFrame(
    processed_h5ad.X,
    schema = processed_h5ad.var_names.to_list()
    ).with_columns([
    pl.Series('cell', processed_h5ad.obs_names.to_list())
])

# move cell column to the first column
h5ad_pl = h5ad_pl.select(['cell'] + h5ad_pl.columns[:-1])

# load in the proportion of significant cells in cell type:
supp1 = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/port/met_scdrs_supp_table/supplementary_table_1.csv')

# get the mean of processed h5ad:
numeric_cols = [col for col, dtype in h5ad_pl.schema.items() if dtype==pl.Float32]
row_mean = h5ad_pl.with_columns(
    pl.mean_horizontal(numeric_cols).alias('row_mean')
)['row_mean']

# load in the MDD result:
mdd_score = pd.read_csv('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.score.gz', sep = '\t')
assert (mdd_score["cell"] == h5ad_pl["cell"].to_pandas()).all()
mdd_score['row_mean'] = row_mean

# find correlation between row_mean and zscore:
mdd_score['norm_score'].corr(mdd_score['row_mean'], method = 'spearman')

# group by cell type and plot:
assert (meta.cell == h5ad_pl["cell"].to_pandas()).all()
meta['row_mean'] = row_mean

# plot supplementary figures:
for cell_class in meta._CellClass.unique():
    cells_in_cell_class = meta[meta._CellClass.isin([cell_class])]
    cells_in_cell_class = cells_in_cell_class.rename(columns={"_MajorType": "MajorType"})
    
    plt.figure(figsize=(10, 6))
    ax = sns.kdeplot(
        data=cells_in_cell_class,
        x="row_mean",
        hue="MajorType",
        common_norm=False,  # Each density integrates to 1 separately
        fill=True,          # Filled density curves
        alpha=0.4           # Transparency
    )
    
    plt.title(f"Density of average mch by cell type: {cell_class}", fontsize=16)
    plt.xlabel("Average mch", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    plt.tight_layout()
    
    # Save to file instead of showing on screen
    output_path = f"/u/home/l/lixinzhe/project-geschwind/plot/{today}_{cell_class}_density_plot_average_mch_by_cell_type.png"
    plt.savefig(output_path, dpi=300)  # high resolution
    plt.close()
