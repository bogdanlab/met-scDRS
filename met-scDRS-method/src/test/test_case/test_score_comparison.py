import pytest
import subprocess
import met_scdrs
from met_scdrs.cli import compare_score
import os

def test_compare_score(tmp_path):
    ###########################################################################################
    ######                                 define tester in                              ######
    ###########################################################################################
    out_folder = tmp_path
    score1 = "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var_length/"
    score2 = "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/v1.1/ges132489_30K_subset/mean_var/"
    plot_path = tmp_path
    
    ###########################################################################################
    ######                                    define process                             ######
    ###########################################################################################
    compare_score(score1, score2, plot_path)
    