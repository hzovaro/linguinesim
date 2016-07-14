############################################################################################
#
# 	File:		mu_min.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Determining the maximum magnitude of stars that can be reliably imaged with
#	LINGUINI.
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
plt.close('all')

# Input parameters


# Seeing
seeing_as = 2	# seeing (arcsec)
sigma_tt_px = 0.5 * seeing_as / SYSTEM_PLATE_SCALE_AS_PX	# tip/tilt standard deviation (px)
crop_tt = np.ceil(3 * sigma_tt_px).astype(int)
detector_size_padded_px = tuple(2*crop_tt+x for x in detector.size_px)

# Image parameters
N_tt = 50			# How many exposures to take for each starfield.
t_exp = 10e-3		# Exposure time
band = 'K'			# Imaging band
n_mags = 11			# Number of star magnitude steps
m = np.linspace(18, 20, n_mags)
coords = np.array([detector.height_px/2 + crop_tt, detector.width_px/2 + crop_tt])	# Star coordinates

# Output parameters
errs = np.zeros((len(m), N_tt))			# Stacking errors
n_errs = np.zeros((len(m)))				# Number of stacking errors for each star magnitude
misaligned_frac = np.zeros((len(m)))	# Fraction of images with stacking errors for each star magnitude

# Make a star field.
for k in range(len(m)):
	print('Simulating star with m = %.2f' % (m[k]))
	starfield_padded_count, starfield_padded_ideal = getStarField(m = m[k], coords = coords, t_exp = t_exp, A_tel = telescope.A_collecting, f_ratio = telescope.f_ratio, l_px_m = detector.l_px_m, detector_size_px = detector_size_padded_px, magnitudeSystem = 'AB', band = band, plotIt = False, tau = telescope.tau, gain = detector.gain, qe = detector.qe)[0:2]

	# Using the count image (computed via Poisson statistics)
	starfield = starfield_padded_count[crop_tt:crop_tt+detector.height_px,crop_tt:crop_tt+detector.width_px]

	# Add tip and tilt.
	starfields_tt, in_idxs = addTurbulence(starfield_padded_count, N_tt, sigma_tt_px, crop_tt)

	# Add noise.
	starfields_noisy = addNoise(starfields_tt, band = band, t_exp = t_exp, plotIt = False)[0]
	etc_output = exposureTimeCalc(band = band, t_exp = t_exp, worstCaseSpider=False, magnitudeSystem = 'AB', surfaceBrightness = m[k])

	# Applying the shift-and-stack routine.
	# Fow now, the reference image is the truth image (review this!)
	starfield_stacked, out_idxs = shiftAndStack(starfields_noisy, image_ref = starfield, plotIt = False)

	# Printing the alignment error.
	n_errs[k], errs[k] = alignmentError(in_idxs, out_idxs)
	# Want to save the fraction of misaligned images 
	misaligned_frac[k] =  n_errs[k] / N_tt

# WANT: a plot of success rate (/100 tries) vs. magnitude
plt.figure(figsize=(FIGSIZE,FIGSIZE))
plt.plot(m, misaligned_frac*100, 'ro')
plt.xlabel('Star AB magnitude')
plt.ylabel('Fraction of misaligned images (\%)')
plt.show()