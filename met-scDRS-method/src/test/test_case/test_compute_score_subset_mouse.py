import pytest
import subprocess
import met_scdrs
from met_scdrs.cli import compute_score
import os

def test_mouse_sim_subset_mean_var_length(tmp_path):
    ###########################################################################################
    ######                                 define tester in                              ######
    ###########################################################################################
    out_folder = tmp_path
    gs_file = '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/small_test_from_publication.gs'
    loaded_gs = met_scdrs.util.load_gs(gs_file)
    
    ###########################################################################################
    ######                                    define process                             ######
    ###########################################################################################
    compute_score(
        h5ad_file = '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad',
        preprocess = True,
        preprocess_method = 'inverse',
        variance_clip = 5,
        h5ad_species = 'mouse',
        cov_file = '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv',
        gs_file = gs_file,
        gs_species = 'human',
        out_folder = out_folder,
        ctrl_match_opt = 'mean_var_length',
        weight_opt = 'inv_std',
        n_ctrl = 1000,
        flag_return_ctrl_raw_score =  True,
        flag_return_ctrl_norm_score = True,
        diagnostic  = True,
        verbose = True
    )
    
    # for each of the result in gs file, make sure the subsequent score output exist:
    for f in loaded_gs.keys():
        score_path = out_folder / f"{f}.score.gz"
        full_score_path = out_folder / f"{f}.full_score.gz"
        assert score_path.exists()
        assert full_score_path.exists()
    