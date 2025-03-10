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
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap

# read in the disease:
trait = sys.argv[1]
aggregation_score_folder = sys.argv[2]
visual_output_dir = sys.argv[3]
common_cortex_axis = sys.argv[4]

# for test:
# trait = "2024-09-16-PASS_MDD_Howard2019-X_MajorType-X_Regionsummary"
# aggregation_score_folder = "/Volumes/misses/research/Bogdan-lab/test/GSE215353-production-aggregate-all/"
# visual_output_dir = "/Volumes/misses/research/Bogdan-lab/spatial-render-out/"
if "true" in sys.argv[4]:
    common_cortex_axis = True
else:
    common_cortex_axis = False

# check if the output exist, if exist, then skip execution:
destination = visual_output_dir + "/movie/" + trait + "/"
if os.path.exists(destination):
    print('interactive object detected, skipping execution!')
else:
    # read in csv for all 75 disease:
    aggregation_scores_folder = aggregation_score_folder
    
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
    output_dir = f"{visual_output_dir}/region-cell-type-aggregated-mean-score/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    for disease, score in aggregates.items():
        score['AIBS_region'] = score['group2.index'].map(AIBS_region_map)
        filtered_score[disease] = score[score['AIBS_region'].isin(AIBS_loading_region['acronym'])]
        filtered_score[disease] = filtered_score[disease].groupby(['group1.index', 'AIBS_region'])['mean_score'].mean().reset_index()
        all_disease_score.extend(filtered_score[disease]['mean_score'].tolist())
        # write the filtered_score out as one of the supplementary table:
        aggregate_out = pandas.DataFrame(filtered_score[disease])
        aggregate_out.to_csv(f"{output_dir}{disease}-region-cell-type-aggregate.csv", index = True)
    
    # compute the normalization factor:
    max_value = max(all_disease_score)
    min_value = min(all_disease_score)
    
    # Create color scale:
    filtered_score = copy.deepcopy(filtered_score)
    for disease, score in filtered_score.items():
        value = filtered_score[disease]['mean_score']
        normalized_values = (value - min_value) / (max_value - min_value)
        filtered_score[disease].loc[:, 'normalized_score'] = normalized_values
    
    # output this color scale into a plot:
    # define green and red color
    green = [0, 1, 0]
    red = [1, 0, 0]
    
    # get a vector of color (unnormalizd by max and min)
    original_color_value = numpy.array(filtered_score[disease]['mean_score'])
    original_color_value = numpy.append(original_color_value, max_value)
    original_color_value = numpy.append(original_color_value, min_value)
    
    # map the color value to the color coding:
    cmap = LinearSegmentedColormap.from_list("GreenRed", [green, red])
    norm = Normalize(vmin=min_value, vmax=max_value)
    fig, ax = plt.subplots()
    
    # Create dummy data and plot (required for color mapping)
    data = numpy.random.rand(10, 10) * 100
    cax = ax.imshow(data, cmap=cmap, norm=norm)  # Use colormap
    plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax = ax, label='Value')
    
    #output the plt to some file: 
    out_dir = visual_output_dir + "interactive/" + trait + "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    plt.savefig(f"{out_dir}{disease}-color-bar.png", dpi = 500)
    
    # subset to the trait score:
    score = filtered_score[trait]
    
    # filter to common cortex regions if the argument is true:
    if common_cortex_axis:
        common_cortex_region = ['V1C', 'V2', 'S1C', 'Pro', 'M1C', 'A5-A7', 'A1C', 'A19', 'A44-A45', 'A46', 'Ig', 'ITG', 'MTG', 'TH-TL', 'A25', 'Idg', 'ACC', 'A38', 'FI', 'Pir', 'MEC', 'LEC']
        common_cortex_mapped = pandas.Series(common_cortex_region).map(AIBS_region_map).tolist()
        # assert set(common_cortex_mapped).issubset(set(score['AIBS_region'])), f"{set(common_cortex_mapped) - set(score['AIBS_region'])} are missing in df['region']"
        score = score[score['AIBS_region'].isin(common_cortex_mapped)]
    
    ### VISUALIZATION #################################################################################
    
    # get cell types:
    cell_types = filtered_score[disease]['group1.index'].unique()
    
    # remove command line output:
    # for each cell type create scene and render:
    for cell_type in cell_types:
        plot_df = score[score['group1.index'] == cell_type]
        unique_regions_count = plot_df['AIBS_region'].nunique()
        print(f"There are {unique_regions_count} unique regions in cell type {cell_type}.")
        
        # Create scene:
        # version 0.1 requried for this to work
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
        out_dir = visual_output_dir + "/movie/" + trait + "/"
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        animation = VideoMaker(scene, out_dir, name = video_name, size='1920x1080')
        animation.make_video(azimuth=3, duration=17, fps=7)
        
        #also output a html interactive object:
        out_dir = visual_output_dir + "interactive/" + trait + "/"
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        interactive_path = out_dir + video_name + ".html"
        scene.export(interactive_path)
        scene.close()
