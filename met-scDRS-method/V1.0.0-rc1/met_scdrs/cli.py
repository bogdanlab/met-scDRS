#!/usr/bin/env python

# import packages:
import fire
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import met_scdrs
from typing import Dict, List
import scanpy as sc
import os
import time

### HEADER ########################################################################################
def get_cli_head():
    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* methylation single-cell disease relevance score (met-scDRS)\n"
    MASTHEAD += "* Version %s\n" % met_scdrs.__version__
    MASTHEAD += "* Xinzhe Li\n"
    MASTHEAD += "* Adapted from scDRS written by Martin Jinye Zhang and Kangcheng Hou\n"
    MASTHEAD += "* UPENN\n"
    MASTHEAD += "* MIT License\n"
    MASTHEAD += "******************************************************************************\n"
    return MASTHEAD

def quote_from_cyberpunk2077():
    # meant for testing script
    print('To This!')

def main():
    fire.Fire()
