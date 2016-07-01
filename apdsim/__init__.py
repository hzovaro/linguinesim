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
#	A package for simulating Lucky Imaging of Nearby Galaxies Undertaken In the Near-Infrared (LINGUINI)
#
#	Copyright (C) 2016 Anna Zovaro
#
###########################################################################################################
#
#	This file is part of lignuini-sim.
#
#	lignuini-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	lignuini-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with lignuini-sim.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################################                                                               
from __future__ import division

# Vega band magnitudes calculated using data from https://www.astro.umd.edu/~ssm/ASTR620/mags.html
VEGA_MAGNITUDE_ZEROPOINT = {
	'J' : 49.46953099,
	'H' : 49.95637318,
	'K' : 50.47441871
}
AB_MAGNITUDE_ZEROPOINT = 48.6
FIGSIZE = 7.5

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

from astropy.io import fits

# Importing detector and telescope properties
import sysparams.cryo as cryo
import sysparams.detector_saphira as detector
import sysparams.telescope_anu23m as telescope

# Importing modules
from apdsim.imutils import *
from apdsim.etcutils import *
from apdsim.etc import *
from apdsim.obssim import *
from apdsim.lisim import *
from apdsim.galsim import *


