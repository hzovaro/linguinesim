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
#	This file is part of linguinesim.
#
#	linguinesim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	linguinesim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with linguinesim.  If not, see <http://www.gnu.org/licenses/>.
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

FILTER_BANDS_M = {
	# [centre wavelength_m, width, min, max]
	# Bands U through I taken from https://en.wikipedia.org/wiki/Photometric_system.
	'U' : [365e-9, 66e-9],
	'B' : [445e-9, 94e-9],
	'V' : [551e-9, 88e-9],
	'R' : [658e-9, 138e-9],
	'I' : [806e-9, 149e-9],
	'J' : [1.250e-6, 0.160e-6, 1.170e-6, 1.330e-6],	# GMTIFS
	'H' : [1.635e-6, 0.290e-6, 1.490e-6, 1.780e-6],	# GMTIFS	
	'K' : [2.200e-6, 0.340e-6, 2.030e-6, 2.370e-6]	# GMTIFS
}
####################################################################################################
" Various packages "
import scipy.constants as constants
import scipy.integrate as integrate
import scipy.optimize as opt
import scipy.signal as signal
import scipy.special as special
from scipy.ndimage.interpolation import shift
import miscutils as mu
import numpy as np
import pdb

# Image processing library
import PIL
from PIL import Image
try:
	from PIL.Image import LANCZOS as RESAMPLE_FILTER
except:
	from PIL.Image import BILINEAR as RESAMPLE_FILTER
	
import os
import sys

# Multithreading/processing packages
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool	# dummy = Threads
from multiprocessing import Pool as ProcPool			# no dummy = Processes

try:
    import pyfftw    
    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(1.0)
    NTHREADS = 2
    print("WARNING: using pyfftw in some routines with NTHREADS = {:d}".format(NTHREADS))
except:
	NTHREADS = 0
	print("WARNING: not using pyfftw; NTHREADS = {:d}".format(NTHREADS))

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary_r')

from astropy.io import fits
from astropy.modeling import models, fitting

# Importing detector and telescope properties
# from linguinesim import sysparams
import linguinesim.sysparams.cryo as cryo
import linguinesim.sysparams.detector_saphira as detector
import linguinesim.sysparams.telescope_anu23m as telescope
import linguinesim.sysparams.sky_sso as sky
SYSTEM_PLATE_SCALE_AS_PX = detector.l_px_m * telescope.plate_scale_as_m
SYSTEM_PLATE_SCALE_RAD_PX = np.deg2rad(SYSTEM_PLATE_SCALE_AS_PX / 3600)
OMEGA_PX_RAD = SYSTEM_PLATE_SCALE_RAD_PX * SYSTEM_PLATE_SCALE_RAD_PX

# Importing modules
from linguinesim.apdsim.imutils import *
from linguinesim.apdsim.etcutils import *
from linguinesim.apdsim.etc import *
from linguinesim.apdsim.obssim import *
from linguinesim.apdsim.lisim import *
from linguinesim.apdsim.galsim import *
from linguinesim.apdsim.starsim import *
from linguinesim.apdsim.satsim import *
from linguinesim.apdsim import fftwconvolve as fftwconvolve



