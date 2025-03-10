### GSE215353-visualization-AIBs.py ###############################################################
# purpose: visualize the mean risk score in each brain disection region for cell type

### DATA LOADING ##################################################################################
# load in libraries:
import pandas
from brainrender import Scene, Animation
from brainrender.video import VideoMaker
import glob
import os
import copy
import re
import gc
import contextlib
import progressbar
import sys

# read in the disease:
trait = sys.argv[1]

# read in csv for all 75 disease:
aggregation_scores_folder = '/Volumes/misses/research/Bogdan-lab/GSE215353-aggregate-zscore/GSE215353-aggregate-all/'

# next for each file inside of the folder, read in the csv file:
aggregates = {}
for csv_path in glob.glob(f'{aggregation_scores_folder}*.csv'):
    # load in each of the disease:
    region = pandas.read_csv(csv_path)
    disease = os.path.basename(csv_path).split('.csv')[0]
    aggregates[disease] = region

### PROCESSING ###################################################################################
# create a region map:
AIBS_region_map = {
    'A19': 'SOG',
    'A1C': 'TTG',
    'A24': 'CgGr',
    'A25': 'SCG',
    'A38': 'TP',
    'A44-A45': 'IFG',
    'A46': 'MFG',
    'A5-A7': 'SPL',
    'ACC': 'ACC',
    'Amy': 'AMY',
    'BL-La': 'BL-La',
    'BM': 'BM',
    'BNST': 'BNST',
    'CA1C-CA2C-CA3C': 'HiT',
    'CA1C-CA2C-CA3C-DGC-CA4C': 'HiT',
    'CA1R-CA2R-CA3R': 'HiB',
    'CA1R-CA2R-CA3R-DGR-CA4R': 'HiB',
    'CBL': 'CBL',
    'CBV': 'CBV',
    'CMN': 'CMN',
    'CaB': 'CaB',
    'Cla': 'Cla',
    'DGC-CA4Cpy': 'HiT',
    'DGR-CA4Rpy': 'HiB',
    'FI': 'FI',
    'GPe': 'GPe',
    'IC': 'IC',
    'IDG': 'IDG',
    'ITG': 'ITG',
    'Idg': 'LIG',
    'Ig': 'SIG',
    'LEC': 'APH',
    'M1C': 'PrCG',
    'MD': 'IC',
    'MEC': 'APH',
    'MTG': 'MTG',
    'NAC': 'NAC',
    'PN': 'Pn',
    'Pir': 'Pir',
    'Pro': 'Pro',
    'Pu': 'Pu',
    'Pul': 'THM',
    'S1C': 'PoCG',
    'SEP': 'SEP',
    'SI': 'SI',
    'Sub': 'HiH',
    'TH-TL': 'PPH',
    'V1C': 'LiG',
    'V2': 'CUN',
    'CEN': 'AMY'
}

# now load in the structure csv from:
AIBS_loading_region = pandas.read_csv('/Volumes/misses/Project/Data/brain-ALLEN-atlas/allen_human_500um_v0.1/structures.csv')

### PROCESSING ####################################################################################
# for each of the disease, we filter to the set of regions that are present in the AIBS 100um atlas:
filtered_score = {}
all_disease_score = []
for disease, score in aggregates.items():
    score['AIBS_region'] = score['group2.index'].map(AIBS_region_map)
    filtered_score[disease] = score[score['AIBS_region'].isin(AIBS_loading_region['acronym'])]
    filtered_score[disease] = filtered_score[disease].groupby(['group1.index', 'AIBS_region'])['mean_score'].mean().reset_index()
    all_disease_score.extend(filtered_score[disease]['mean_score'].tolist())

# compute the normalization factor:
max_value = max(all_disease_score)
min_value = min(all_disease_score)

# Create color scale:
filtered_score = copy.deepcopy(filtered_score)
for disease, score in filtered_score.items():
    value = filtered_score[disease]['mean_score']
    normalized_values = (value - min_value) / (max_value - min_value)
    filtered_score[disease].loc[:, 'normalized_score'] = normalized_values

# subset to the trait score:
score = filtered_score[trait]

### VISUALIZATION #################################################################################
# define color:
green = [0, 1, 0]
red = [1, 0, 0]

# get cell types:
cell_types = filtered_score[disease]['group1.index'].unique()

# remove command line output:
# for each cell type create scene and render:
for cell_type in cell_types:
    plot_df = score[score['group1.index'] == cell_type]
    
    # Create scene:
    scene = Scene(title = cell_type, atlas_name = "allen_human_500um")
    
    # for each of the region, add the normalized score into the rendering object:
    for index, entry in plot_df.iterrows():
        region = entry['AIBS_region']
        value = entry['normalized_score']
        color = [green[i] + (red[i] - green[i]) * value for i in range(3)]
        scene.add_brain_region(region, alpha=0.5, color=color)
    
    # Create video render:
    video_name = trait + '-' + cell_type
    video_name = re.sub(r'/', '-', video_name)
    out_dir = "/Volumes/misses/research/Bogdan-lab/GSE215353-aggregate-zscore/GSE215353-aggregate-all-renders/movie/" + trait + "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    animation = VideoMaker(scene, out_dir, name = video_name, size='1920x1080')
    animation.make_video(azimuth=3, duration=17, fps=7)
    
    #also output a html interactive object:
    out_dir = "/Volumes/misses/research/Bogdan-lab/GSE215353-aggregate-zscore/GSE215353-aggregate-all-renders/interactive/" + trait + "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    interactive_path = out_dir + video_name + ".html"
    scene.export(interactive_path)
    scene.close()
