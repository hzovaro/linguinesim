###########################################################################################################
#   ________  ________  ________  ________  ___  _____ ______      
#  |\   __  \|\   __  \|\   ___ \|\   ____\|\  \|\   _ \  _   \    
#  \ \  \|\  \ \  \|\  \ \  \_|\ \ \  \___|\ \  \ \  \\\__\ \  \   
#   \ \   __  \ \   ____\ \  \ \\ \ \_____  \ \  \ \  \\|__| \  \  
#    \ \  \ \  \ \  \___|\ \  \_\\ \|____|\  \ \  \ \  \    \ \  \ 
#     \ \__\ \__\ \__\    \ \______/|____\_\  \ \__\ \__\    \ \__\
#      \|__|\|__|\|__|     \|______/\__________\|__|\|__|     \|__|
#                                   \|_________|                   
#
#	File:		__init__.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A package for simulating Lucky Imaging of Nearby Galaxies (Usefully) In the Near-Infrared (LINGUINI)
#
############################################################################################################                                                               
from __future__ import division

# Required packages
import scipy.constants as constants
import scipy.integrate as integrate
import scipy.optimize as opt
import scipy.signal as signal
import scipy.special as special
from scipy.ndimage.interpolation import shift

import numpy as np

import pdb

import PIL
from PIL import Image

import os

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary')
figsize = 7.5

from astropy.io import fits

# Importing detector and telescope properties
import sysparams.cryo as cryo
import sysparams.detector_saphira as detector
import sysparams.telescope_anu23m as telescope

# Importing modules
from apdsim.imutils import *
from apdsim.etc import *
from apdsim.obssim import *
from apdsim.lisim import *

