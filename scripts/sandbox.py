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

wavelength = 2e-6
R = 50*25
D = 2.5
A_tel = np.pi * D * D / 4
f_ratio = R / D
l_px_m = 20e-6
width = 500
height = 500
A_detector = height * width * l_px_m * l_px_m
detector_size_px = (height,width)

# psf, P_0, I_0 = psfKernel(wavelength, f_ratio, l_px_m, detector_size_px,
# 	plotIt=True)
# psf, P_0, I_0 = psfKernel(2.2e-6, telescope.f_ratio, detector.l_px_m, detector.size_px,
# 	plotIt=True)
# print 'Sum = ',sum(psf.flatten()) * detector.l_px_m * detector.l_px_m
# print 'P_0 = ',P_0
# print 'I_0 = ',I_0

# star = getStar(funtype, coords, wavelength, f_ratio, l_px_m, detector_size_px,
# 	plotIt=False)

image_count, starfield, m, coords = getStarField(N_stars = 5, A_tel = telescope.A_collecting, f_ratio = telescope.f_ratio, l_px_m = detector.l_px_m, detector_size_px = detector.size_px, magnitudeSystem = 'AB', band = 'K', plotIt = True)
