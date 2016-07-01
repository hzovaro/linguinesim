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
#
#	TO DO:
#	- Check consistency of np.log usage!
#	- better way of calculating b_n
#	- weird pixelation effect: occurs during call to image.rotate()
#	- Ultimate goal: what is the minimum surface brightness of galaxies that we can reliably
#	image with LINGUINI? 
#	- Image directly onto the APD 
#
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
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
z = 0.095		# Redshift
D_A_Mpc = distances(z)['D_A_Mpc']	# Angular diameter distance of galaxy (Mpc)
D_A_kpc = D_A_Mpc * 1e3 			# Angular diameter distance of galaxy (kpc)
R_e = 4			# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
R_trunc = 50	# Truncation radius of disc
mu_e = 15		# Surface brightness at the half-light radius (AB mag/arcsec^2)
n = 4			# Sersic index
i_deg = 45		# inclination (degrees)
theta_deg = 0	# rotation (degrees)

# Observing parameters:
t_exp = 100e-3	# exposure time (s)
band = 'K'		# imaging band 

# 1. Simulating a galaxy.
# In units of kpc...

# Distance subtended by one pixel a distance D_A_kpc away given our plate scale.
dR_kpc = np.deg2rad(SYSTEM_PLATE_SCALE_AS_PX / 3600) * D_A_kpc
R_max = dR_kpc * detector.width_px / 2

R, dR, I_map, mu_map = sersic2D(n=n, R_e=R_e, R_trunc=R_trunc, mu_e=mu_e,
	R_max = R_max,
	gridsize = detector.width_px,
	i_rad=np.deg2rad(i_deg), 
	theta_rad=np.deg2rad(theta_deg),
	plotIt=True)

