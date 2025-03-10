### contex-common-axis-visualization.py ###########################################################
# purpose: visualize the common axis in the cortex

### PREAMBLE ######################################################################################
# load in libraries
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

### VISUALIZATION #################################################################################
# Create a scene
scene = Scene(title = 'common cortex axis', atlas_name = "allen_human_500um")

# add Area 44-45 and area 46
scene.add_brain_region('IFG', alpha = 1, color = '#8da0cb') # inferior frontal gyrus - BA44-45 - prefrontal cortex
scene.add_brain_region('MFG', alpha = 1, color = '#fc8d62') # middle frontal gyrus - BA46
scene.add_brain_region('ITG', alpha = 1, color = '#66c2a5') # Inferior temproal gyrus
scene.add_brain_region('MTG', alpha = 1, color = '#e78ac3') #middle temporal gyrus
scene.add_brain_region('LiG', alpha = 1, color = '#fdc086') # lingual gyrus - visual cortex

# render scene:
scene.render(camera = 'sagittal', interactive = False)
scene.screenshot(name="IFG-MFG-ITG-MTG-LiG-render") 
scene.close()
