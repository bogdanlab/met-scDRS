### gse132489_normality_check.py ##################################################################
# purpose: check the normality assumption for the gse132489 30K run before and after normalization

### PREAMBLE ######################################################################################
# load in packages:
import pandas as pd
import re
import os
from tqdm import tqdm
from collections import defaultdict
import numpy as np

# load in all sets of null distributions:
score_paths = [
    "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/arcsine/null_distribution/",
    "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/library/null_distribution/",
    "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/logit/null_distribution/",
    "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/untransformed/null_distribution/"
]

# load in the cell meta:
meta_data = pd.read_csv('/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv', sep = ',', index_col = 0)

# identify the files that are of metric from the list of dirs:
normalization_metrics = defaultdict(dict)
for path in score_paths:
    # get compute mode:
    mode = re.sub("/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/", '', path)
    mode = re.sub("/null_distribution/", "", mode)
    
    # identify files to read in 
    metric_files = [file for file in os.listdir(path) if file.endswith('_metric.txt')]
    significant_count = []
    marginal_lrt = []
    
    # for each of the disease, read in the summary:
    for metric_file in tqdm(metric_files):
        disease = re.sub('_calibration_p_metric.txt', '', metric_file)
        metric_path = os.path.join(path, metric_file)
        metric = pd.read_csv(metric_path, sep = '\t', index_col = 0)
        normalization_metrics[mode][disease] = metric
        
        # calculate the number of significant nulls that are not normal:
        significant_count.append(np.sum(metric.ctrl_raw_adjusted < 0.05))
    
    # print the average number of significant count across traits:
    print(f'average number of significant nulls after correction across traits {np.mean(significant_count)}')
