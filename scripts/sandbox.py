############################################################################################
#
# 	File:		sandbox.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A 'sandbox' for testing various features. Play safe kids!
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
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
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
plt.close('all')

# Input parameters
N_stars = 5
N_tt = 5
t_exp = 10e-3
band = 'K'
sigma_tt_px = 20
crop_tt = 50
m_max = 20
m_min = 10

# Make a star field.
image_count, starfield_padded, m, coords = getStarField(N_stars = N_stars, m_min = m_min, m_max = m_max, A_tel = telescope.A_collecting, f_ratio = telescope.f_ratio, l_px_m = detector.l_px_m, detector_size_px = tuple(2*crop_tt+x for x in detector.size_px), magnitudeSystem = 'AB', band = band, plotIt = True)
starfield = starfield_padded[crop_tt:crop_tt+detector.height_px,crop_tt:crop_tt+detector.width_px]

# Add tip and tilt.
starfields_tt, in_idxs = addTurbulence(starfield_padded, N_tt, sigma_tt_px, crop_tt)

# Add noise.
starfields_noisy, etc_output = addNoise(starfields_tt, band = band, t_exp = t_exp)

# Applying the shift-and-stack routine.
starfield_stacked, out_idxs = shiftAndStack(starfields_tt, image_ref = starfield, plotIt = True)

# Printing the alignment error.
printAlignmentError(in_idxs, out_idxs)