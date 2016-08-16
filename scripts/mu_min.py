############################################################################################
#
# 	File:		mu_min.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Determining the minimum surface brightness of galaxies that can be reliably imaged in
#	LINGUINI.
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
#
#	This file is part of lingiune-sim.
#
#	lingiune-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	lingiune-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with lingiune-sim.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
plt.close('all')
"""
	STEPS:
	1. Simulate a galaxy. For now, make the effective surface brightness mu_e the only variable parameter. For the morphology assume a combined disc (n = 1) + bulge (n = 4) profile. For the effective radius assume the mean value for SAMI galaxies (R_e ~ 4-5 kpc). Assume an appropriate exposure time (matched to tau_0) for the given imaging band.
	2. Resize the image to the detector.
	3. Convolve with the PSF of the 2.3 m.
	4. Make N_tt copies of the image with randomised tip and tilt.
	5. Add noise to each copy.
	6. Shift-and-stack.

"""

" Inputs "
# Galaxy:
z = 0.01		# Redshift
D_A_Mpc = distances(z)['D_A_Mpc']	# Angular diameter distance of galaxy (Mpc)
D_A_kpc = D_A_Mpc * 1e3 			# Angular diameter distance of galaxy (kpc)
R_e = 4			# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
R_trunc = np.inf	# Truncation radius of disc
mu_e = 13		# Surface brightness at the half-light radius (AB mag/arcsec^2)
n = 4			# Sersic index
i_deg = 75		# inclination (degrees)
theta_deg = 0	# rotation (degrees)

# Observing parameters:
t_exp = 100e-3	# exposure time (s)
band = 'K'		# imaging band 
N_exp = 50		# number of exposures
seeing_as = 2	# seeing (arcsec)
sigma_tt_px = 0.5 * seeing_as / SYSTEM_PLATE_SCALE_AS_PX	# tip/tilt standard deviation (px)

# 1, 2. Simulating a galaxy.
# Distance subtended by one pixel a distance D_A_kpc away given our plate scale.
dR_kpc = np.deg2rad(SYSTEM_PLATE_SCALE_AS_PX / 3600) * D_A_kpc
R_max = dR_kpc * detector.width_px / 2

R, dR, F, mu = sersic2D(n=n, R_e=R_e, R_trunc=R_trunc, mu_e=mu_e,
	R_max = R_max,
	gridsize = detector.width_px,
	i_rad=np.deg2rad(i_deg), 
	theta_rad=np.deg2rad(theta_deg),
	wavelength_m = FILTER_BANDS_M[band][0],
	zeropoint = AB_MAGNITUDE_ZEROPOINT,
	plotIt=True)

# 3. Convolving the image with the PSF of the 2.3 m.
# 3a. Convert to counts.
countrate_truth = flux2countRate(F = F, A_tel = telescope.A_collecting, plate_scale_as_px = SYSTEM_PLATE_SCALE_AS_PX, band = band, magnitudeSystem = 'AB', tau = telescope.tau, qe = detector.qe, gain = detector.gain)
image_truth = expectedCount2count(countrate_truth, t_exp)

# 3b. Diffraction limit the image
countrate_difflim = getDiffractionLimitedImage(image_truth = countrate_truth, wavelength = FILTER_BANDS_M[band][0], f_ratio = telescope.f_ratio, l_px_m = detector.l_px_m, detector_size_px = detector.size_px, plotIt = True)
image_difflim = expectedCount2count(countrate_difflim, t_exp)

# 4. Making N_exp copies with randomised tip and tilt.
images_tt, tt_idxs = addTipTilt(image_difflim, N_tt = N_exp, sigma_tt_px = sigma_tt_px)

# 5. Adding noise to each copy.
images_tt_noisy, etc_output = addNoise(images = images_tt, band = band, t_exp = t_exp, plotIt = True)

# 6. Shifting-and-stacking.
image_stacked = xcorrShiftAndStack(images = images_tt_noisy, image_ref = image_difflim, plotIt = True)

exportGalaxyFITSFile(image_stacked, n, R_e, mu_e, z, R_trunc, i_deg, band, seeing_as, t_exp, N_exp, relpath='../galfit')