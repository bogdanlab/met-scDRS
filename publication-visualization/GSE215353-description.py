### GSE215353-description.R #######################################################################
# purpose: describe the GSE215353 dataset in how the MCH distribute by cell type, regions ...

### PREAMBLE ######################################################################################
# load in libraries: 
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime

# load in the dataset:
data_path = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mch.h5ad'
adata = sc.read_h5ad(data_path)
meta = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv')
system_date = datetime.date.today()

### PROCESS #######################################################################################
# compute the row averge:
row_average = np.mean(adata.X, axis = 1)
column_average = np.mean(adata.X, axis = 0)

# check if the cell name in the adata is identical to meta
assertion = adata.obs_names == meta.cell
assert assertion.all()

# since the cell names are identical in order, we can simply merge the average of methylation into meta:
meta_agg = pd.concat([meta, pd.DataFrame(row_average)], axis = 1)
meta_agg.columns = meta.columns.to_list() + ['avg_mch']

# make plot on the density:
g = sns.FacetGrid(meta_agg, hue = '_MajorType', col='_CellClass', height=5, aspect=1, sharex=True, sharey=True)
g.map(sns.kdeplot, 'avg_mch', fill=False, alpha=0.6)

# Show the plot
plt.savefig(f'/u/home/l/lixinzhe/project-geschwind/plot/{system_date}_density_plot.png', dpi = 400)
plt.title("Density Plot")
plt.legend(title="Cell Type")
plt.xlabel("Value")
plt.ylabel("Density")
plt.show()
plt.close('all')

# for each of the cell class, plot out the the cell type:
exc_cell_type = meta[meta._CellClass == 'Excitatory Neurons']['_MajorType'].unique()
inh_cell_type = meta[meta._CellClass == 'Inhibitory and Subcortical Neurons']['_MajorType'].unique()
nn_cell_type = meta[meta._CellClass == 'Non-neuronal Cells']['_MajorType'].unique()
for num,cell_type in enumerate(exc_cell_type):
    subset = meta_agg[meta_agg['_MajorType'] == cell_type]
    kde = sns.kdeplot(data=subset, x='avg_mch', label=None, alpha=0.6)
    
    # get coordinates for annotation:
    x, y = kde.lines[num].get_data()
    max_idx = np.argmax(y)
    plt.annotate(
        cell_type,  # Text label
        xy=(x[max_idx], y[max_idx]),  # Location of the peak
        fontsize=10,
        color='black'
    )

plt.title('Excitatory neurons mch distribution')
plt.xlabel("average mch level")
plt.xlim(0.8, 1)
plt.ylabel("Density")
plt.savefig(f'/u/home/l/lixinzhe/project-geschwind/plot/{system_date}_exc_density_plot_arrow.png', dpi = 400)
plt.close('all')

# Next, look at the plot for inhibitory neurons:
for num,cell_type in enumerate(inh_cell_type):
    subset = meta_agg[meta_agg['_MajorType'] == cell_type]
    kde = sns.kdeplot(data=subset, x='avg_mch', label=None, alpha=0.6)
    
    # get coordinates for annotation:
    x, y = kde.lines[num].get_data()
    max_idx = np.argmax(y)
    plt.annotate(
        cell_type,  # Text label
        xy=(x[max_idx], y[max_idx]),  # Location of the peak
        fontsize=10,
        color='black'
    )

plt.title('Inh neurons mch distribution')
plt.xlim(0.8, 1)
plt.xlabel("average mch level")
plt.ylabel("Density")
plt.savefig(f'/u/home/l/lixinzhe/project-geschwind/plot/{system_date}_inh_density_plot_arrow.png', dpi = 400)
plt.close('all')

# finally, look at the plot for non neuronal cells:
for num,cell_type in enumerate(nn_cell_type):
    subset = meta_agg[meta_agg['_MajorType'] == cell_type]
    kde = sns.kdeplot(data=subset, x='avg_mch', label=None, alpha=0.6)
    
    # get coordinates for annotation:
    x, y = kde.lines[num].get_data()
    max_idx = np.argmax(y)
    plt.annotate(
        cell_type,  # Text label
        xy=(x[max_idx], y[max_idx]),  # Location of the peak
        fontsize=10,
        color='black'
    )

plt.title('non-neurons mch distribution')
plt.xlim(0.8, 1)
plt.xlabel("average mch level")
plt.ylabel("Density")
plt.savefig(f'/u/home/l/lixinzhe/project-geschwind/plot/{system_date}_non_neuronal_density_plot_arrow.png', dpi = 400)
plt.close('all')
