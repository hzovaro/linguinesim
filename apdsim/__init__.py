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
from __future__ import division, print_function

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

# Importing modules
from linguinesim.apdsim.constants import *
from linguinesim.apdsim.imutils import *
from linguinesim.apdsim.etcutils import *
from linguinesim.apdsim.etc import *
from linguinesim.apdsim.obssim import *
from linguinesim.apdsim.ossim import *
from linguinesim.apdsim.lisim import *
from linguinesim.apdsim.galsim import *
from linguinesim.apdsim.starsim import *
from linguinesim.apdsim.satsim import *
from linguinesim.apdsim import fftwconvolve as fftwconvolve



