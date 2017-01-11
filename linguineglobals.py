	################################################################################
#
# 	File:		linguineglobals.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Some useful constants.
#
#	Copyright (C) 2016 Anna Zovaro
#
################################################################################
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
################################################################################
from __future__ import division, print_function
import astropy.constants

################################################################################
# Vega band magnitudes calculated using data from 
# https://www.astro.umd.edu/~ssm/ASTR620/mags.html
# See also http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
VEGA_MAGNITUDE_ZEROPOINT = {
	'J' : 49.46953099,
	'H' : 49.95637318,
	'K' : 50.47441871
}
AB_MAGNITUDE_ZEROPOINT = 48.6

FILTER_BANDS_M = {
	# [centre wavelength_m, width, min, max]
	# Bands U through I taken from https://en.wikipedia.org/wiki/Photometric_system.
	'U' : [365e-9, 66e-9, 0, 0],
	'B' : [445e-9, 94e-9, 0, 0],
	'V' : [551e-9, 88e-9, 0, 0],
	'R' : [658e-9, 138e-9, 0, 0],
	'I' : [806e-9, 149e-9, 0, 0],
	'J' : [1.250e-6, 0.160e-6, 0, 0],	# GMTIFS
	'H' : [1.635e-6, 0.290e-6, 0, 0],	# GMTIFS	
	'K' : [2.200e-6, 0.340e-6, 0, 0]	# GMTIFS
}
# Calculating filter endpoints
for key in FILTER_BANDS_M:
	FILTER_BANDS_M[key][2] = FILTER_BANDS_M[key][0] - 0.5 * FILTER_BANDS_M[key][1]
	FILTER_BANDS_M[key][3] = FILTER_BANDS_M[key][0] + 0.5 * FILTER_BANDS_M[key][1]

AL_REFLECTIVITY = 0.909
AL_EMISSIVITY = 1 - AL_REFLECTIVITY

# TENTH_AIRY_RING = 10.25		# multiples of N_os corresponding to the 10th Airy ring in a diffraction-limited PSF

# Solar properties
T_SUN_K = 5777 								# Temperature (K)
R_SUN_M = astropy.constants.R_sun.value 	# Radius (m)
DIST_SUN_M = astropy.constants.au.value		# 1 AU (Distance from Earth's centre (m))