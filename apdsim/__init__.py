####################################################################################################
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
####################################################################################################
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
####################################################################################################                                                             
from __future__ import division

####################################################################################################
" Global variables "
# Vega band magnitudes calculated using data from https://www.astro.umd.edu/~ssm/ASTR620/mags.html
VEGA_MAGNITUDE_ZEROPOINT = {
	'J' : 49.46953099,
	'H' : 49.95637318,
	'K' : 50.47441871
}
AB_MAGNITUDE_ZEROPOINT = 48.6

# Plotting
FIGSIZE = 5
COLORBAR_FRACTION = 0.046
COLORBAR_PAD = 0.04

# Near-IR filter bands
FILTER_BANDS_UM = {
	# [centre wavelength_m, width, min, max]
	'J' : [1.250, 0.160, 1.170, 1.330],
	'H' : [1.635, 0.290, 1.490, 1.780],
	'K' : [2.200, 0.340, 2.030, 2.370]
}
FILTER_BANDS_M = {
	# [centre wavelength_m, width, min, max]
	'J' : [1.250e-6, 0.160e-6, 1.170e-6, 1.330e-6],
	'H' : [1.635e-6, 0.290e-6, 1.490e-6, 1.780e-6],
	'K' : [2.200e-6, 0.340e-6, 2.030e-6, 2.370e-6]
}
####################################################################################################
" Various packages "
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

import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary')

from astropy.io import fits
from astropy.modeling import models, fitting

# Importing detector and telescope properties
import sysparams.cryo as cryo
import sysparams.detector_saphira as detector
import sysparams.telescope_anu23m as telescope
import sysparams.sky_sso as sky
SYSTEM_PLATE_SCALE_AS_PX = detector.l_px_m * telescope.plate_scale_as_m
SYSTEM_PLATE_SCALE_RAD_PX = np.deg2rad(SYSTEM_PLATE_SCALE_AS_PX / 3600)
OMEGA_PX_RAD = SYSTEM_PLATE_SCALE_RAD_PX * SYSTEM_PLATE_SCALE_RAD_PX

# Importing modules
from apdsim.imutils import *
from apdsim.etcutils import *
from apdsim.etc import *
from apdsim.obssim import *
from apdsim.lisim import *
from apdsim.galsim import *
from apdsim.starsim import *


