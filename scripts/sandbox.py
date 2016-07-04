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
#
#	TO DO:
#	- Check consistency of np.log usage!
#	- better way of calculating b_n
#	- Ultimate goal: what is the minimum surface brightness of galaxies that we can reliably
#	image with LINGUINI?
#
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
# plt.close('all')

# Simulating an image of a galaxy with a Sersic profile.
# Disc: truncated exponential.
# Bulge: de Vaucoeleurs profile. 

" Inputs "
# Galaxy:
R_e = 4	# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
R_trunc = 50	# Truncation radius of disc
mu_e = 15	# Surface brightness at the half-light radius (AB mag/arcsec^2)
n = 4		# Sersic index
z = 0.095	# Redshift
i_deg = 45		# inclination (degrees)
theta_deg = 0	# rotation (degrees)
alpha = 30	# rotation (degrees)
band = 'K'

# In units of kpc...
R, dR, I_map, mu_map = sersic2D(n=n, R_e=R_e, R_trunc=R_trunc, mu_e=mu_e, plotIt=True, 
	gridsize=1000,
	i_rad=np.deg2rad(i_deg), 
	theta_rad=np.deg2rad(theta_deg))

# Converting to units of pixels on our detector...
D_A_Mpc = distances(z)['D_A_Mpc']
# Angle (in arcsec) that each pixel on our imaginary detector subtends on the sky at redshift z.
plate_scale_as = np.rad2deg(2 * np.arctan(dR / 2 / (D_A_Mpc * 1e3))) * 3600

# Resizing to detector AFTER converting to electron count
flux_map = getElectronCount(mu = mu_map, 
	band = band,
	A_tel = telescope.A_collecting,
	plate_scale_as_px = plate_scale_as,
	tau = telescope.tau,
	qe = detector.qe,
	gain = detector.gain,
	magnitudeSystem = 'AB'
	)

# Resizing to detector BEFORE converting to electron count
mu_map_apd = resizeImagesToDetector(images_raw = mu_map, 
	source_plate_scale_as = plate_scale_as, 
	dest_detector_size_px = (detector.height_px, detector.width_px),
	dest_plate_scale_as = SYSTEM_PLATE_SCALE_AS_PX)

flux_map_apd_before = getElectronCount(mu = mu_map_apd, 
	band = band,
	A_tel = telescope.A_collecting,
	plate_scale_as_px = SYSTEM_PLATE_SCALE_AS_PX,
	tau = telescope.tau,
	qe = detector.qe,
	gain = detector.gain,
	magnitudeSystem = 'AB'
	)

extent_imag_detector_as = [-plate_scale_as*mu_map.shape[1]/2,plate_scale_as*mu_map.shape[1]/2,-plate_scale_as*mu_map.shape[0]/2,plate_scale_as*mu_map.shape[0]/2]
extent_apd_detector_as = [-SYSTEM_PLATE_SCALE_AS_PX*mu_map_apd.shape[1]/2,SYSTEM_PLATE_SCALE_AS_PX*mu_map_apd.shape[1]/2,-SYSTEM_PLATE_SCALE_AS_PX*mu_map_apd.shape[0]/2,SYSTEM_PLATE_SCALE_AS_PX*mu_map_apd.shape[0]/2]

# PLOTTING FLUX AND MU MAPS ON EACH DETECTOR
plt.figure(figsize=(2*FIGSIZE, 2*FIGSIZE))
plt.subplot(2,2,1)
plt.imshow(mu_map, extent = extent_imag_detector_as)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title(r'$\mu(R)$ (imaginary detector)')
plt.xlabel(r'arcsec')
plt.ylabel(r'arcsec')
plt.subplot(2,2,2)
plt.imshow(mu_map_apd, extent = extent_apd_detector_as)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title(r'$\mu(R)$ (eAPD detector)')
plt.xlabel(r'arcsec')
plt.ylabel(r'arcsec')
plt.subplot(2,2,3)
plt.imshow(flux_map, extent = extent_imag_detector_as)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title(r'Electron count (imaginary detector)')
plt.xlabel(r'arcsec')
plt.ylabel(r'arcsec')
plt.subplot(2,2,4)
plt.imshow(flux_map_apd_before, extent = extent_apd_detector_as)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title(r'Electron count (eAPD detector)')
plt.xlabel(r'arcsec')
plt.ylabel(r'arcsec')
plt.show()

# exportFITSFile(image_in_array = flux_map_apd_before, fname = 'flux_data', clobber = True, headerData = {'EXPTIME' : 1.0})

