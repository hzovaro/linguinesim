############################################################################################
#
# 	File:		starsim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	For simulating images of stars. 
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

# Inputs:
#	- star magnitude & magnitude system
#	- star coordinates (in arcsec relative to corner of image)
# Outputs:
#	- an image of the star as formed by the telescope and detector pair.

def getStarField(m, coords, A_tel, f_ratio, l_px_m, detector_size_px,
	tau = 1,
	qe = 1,
	gain = 1,
	magnitudeSystem = None,
	wavelength_m = None, 
	bandwidth_m = None, 
	band = None, 
	plotIt=False):
	
	# Figuring out the value of the pixel:
	# If we were imaging onto our detector with an infinitely large telescope, then all of the light falling on the telescope aperture would be imaged onto a single pixel.
	# So, we need to get the flux (in units of electrons per pixel) given the magnitude of the star.
	# Then, we make a normalised PSF grid (s.t. P = 1) with the appropriate plate scale, and simply multiply this by the flux. 
	# THEN we convert to an expected count.
	N_stars = len(m)
	if len(m) != coords.shape[0]:
		print 'ERROR: the length of the star magnitudes vector must be equal to that of the first dimension of the coords array!'
		return

	# PSF.
	detector_height_px, detector_width_px = detector_size_px[0:2]
	if wavelength_m == None and band != None:
		wavelength_m = FILTER_BANDS_M[band][0]
		bandwidth_m = FILTER_BANDS_M[band][1]
	psf = pointSpreadFunction(wavelength_m, f_ratio, l_px_m, tuple(2*x for x in detector_size_px), plotIt=True)
	# 'Truth' image (telescope with infinite resolution)
	image_truth = np.zeros(detector_size_px)

	for k in range(N_stars):
		# Total expected electron count from the star.
		Sigma_electrons = surfaceBrightness2countRate(mu = m[k], A_tel = A_tel, tau = tau, qe = qe, gain = gain, magnitudeSystem = magnitudeSystem, wavelength_m = wavelength_m, bandwidth_m = bandwidth_m)		
		image_truth[coords[k,0], coords[k,1]] = Sigma_electrons

	# Convolving the truth image with the PSF.
	image_expected = signal.fftconvolve(psf, image_truth, mode='same')[detector_height_px//2 + detector_height_px%2:detector_height_px+detector_height_px//2 + detector_height_px%2 - 1, detector_width_px//2 + detector_width_px%2:detector_width_px+detector_width_px//2 + detector_width_px%2 - 1]
	if min(image_expected.flatten()) < 0:
		print 'WARNING: the result of the convolution contains negative pixels!'
		image_expected -= min(image_expected.flatten())
	image_count = expectedCount2count(image_expected)

	if plotIt:
		plt.figure(figsize=(3*FIGSIZE,FIGSIZE))
		plt.suptitle('Starfield')
		plt.subplot(1,3,1)
		plt.imshow(image_truth, norm=LogNorm())
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Truth image')
		plt.subplot(1,3,2)
		plt.imshow(image_expected, norm=LogNorm())
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Expected electron flux')
		plt.subplot(1,3,3)
		plt.imshow(image_count, norm=LogNorm())
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Count')
		plt.show()

	return image_count, image_expected, image_truth



