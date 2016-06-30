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
#	This file is part of apd-sim.
#
#	apd-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	apd-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with apd-sim.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################
#
#	TO DO:
#	- Check consistency of np.log usage!
#	- Conversion to/from I_e using mu_e: units?
#	- better way of calculating b_n
#
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
plt.close('all')

# Simulating an image of a galaxy with a Sersic profile.
# Disc: truncated exponential.
# Bulge: de Vaucoeleurs profile. 

" Inputs "
# Galaxy:
R_e = 4	# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
R_trunc = np.inf	# Truncation radius of disc
mu_e = 12	# Surface brightness at the half-light radius (AB mag/arcsec^2)
n = 1
z = 0.095	# Redshift
i = 45		# inclination (degrees)
alpha = 30	# rotation (degrees)
band = 'K'

# In units of kpc...
R, dR, I_map, mu_map = sersic2D(n=n, R_e=R_e, R_trunc=R_trunc, mu_e=mu_e, plotIt=True, 
	i_rad=np.deg2rad(30), 
	theta_rad=np.deg2rad(20))

# Converting to units of pixels on our detector...
# First, need to calculate the mapping from parsecs-->arcsec.
D_A_Mpc = distances(z)['D_A_Mpc']
plate_scale_as = np.rad2deg(2 * np.arctan(dR / 2 / (D_A_Mpc * 1e3))) * 3600

flux_map = getPhotonFlux(surfaceBrightness=mu_map, 
	wavelength_eff=telescope.filter_bands_m[band][0], 
	bandwidth=telescope.filter_bands_m[band][1], 
	plate_scale_as=telescope.plate_scale_as_m*detector.l_px_m, 
	A_tel=telescope.A_collecting)

flux_map_resized = resizeImagesToDetector(flux_map, 
	plate_scale_as, 
	(detector.height_px, detector.width_px),
	detector.l_px_m, 
	detector.l_px_m * telescope.plate_scale_as_m)

plt.figure(figsize=(2*figsize,figsize))
plt.subplot(1,2,1)
plt.imshow(flux_map)
plt.colorbar()
plt.title('2D flux map (imaginary detector, plate scale = %.2f arcsec/pixel)' % plate_scale_as)
plt.subplot(1,2,2)
plt.imshow(flux_map_resized)
plt.colorbar()
plt.title('2D flux map (imaginary detector, plate scale = %.2f arcsec/pixel)' % plate_scale_as)
plt.show()

# So the galaxy has a certain redshift and a certain half-light radius
# We want our process to be: Light profile in kpc --> Light profile in arcsec --> Light profile in pixels.
# So make the Sersic profile in units of linear distance (say pc or kpc), then use the redshift to convert 1 kpc into arcseconds acros the sky (use cosmo_calc.py for this), then rescale the Sersic profile to arcseconds, then rescale to pixels on the detector.
# plt.figure(figsize=(2*figsize, figsize))
# for n in range(1,5):
# 	R, mu, I = sersic(n=n, R_e=10, mu_e=mu_e, R=np.linspace(0, 10*R_e, 100))
# 	plt.subplot(1,2,1)
# 	plt.plot(R/R_e, mu, color = (0, n/5, 1), label = r'$n=%d$' % n)
# 	plt.subplot(1,2,2)
# 	plt.plot(R/R_e, I, color = (0, n/5, 1), label = r'$n=%d$' % n)

# plt.subplot(1,2,1)
# plt.plot([1,1], [min(mu), max(mu)], 'k--')
# plt.title(r'Surface brightness vs. radius, $n = %d$' % n)
# plt.xlabel(r'Radius $R/R_e$')
# plt.ylabel(r'Surface brightness $\mu(R)$ (magnitudes/arcsec$^2$)')
# plt.legend()
# plt.gca().invert_yaxis()

# plt.subplot(1,2,2)
# plt.plot([1,1], [0, max(I)], 'k--')
# plt.title(r'Intensity vs. radius, $n = %d$' % n)
# plt.xlabel(r'Radius $R/R_e$')
# plt.ylabel(r'Intensity $I(R)$')
# plt.yscale('log')
# plt.legend()
# plt.show()


