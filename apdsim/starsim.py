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

def getStarField(A_tel, f_ratio, l_px_m, detector_size_px,
	t_exp = 1,			# Exposure time (s).
	N_stars = None,			# Number of random stars to generate.
	m = None, 				# Vector of star magnitudes.
	coords = None, 			# Coordinates of stars.
	m_min = 20,			# Minimum star magnitude.
	m_max = 5,			# Maximum star magnitude.
	tau = 1,
	qe = 1,
	gain = 1,
	magnitudeSystem = None,
	wavelength_m = None, 
	bandwidth_m = None, 
	band = None, 	
	verbose=False,		# If True, prints a table showing the respective magnitudes and coordinates of each star in the field.
	detectorSaturation = np.inf,
	plotIt=False):
	"""
		Returns a simulated image of a star field imaged through an optical system at a given wavelength_m or in an imaging band with a specified centre wavelength_m and bandwidth with a given collecting area, f ratio, pixel size and detector dimensions. The throughput, QE and gain of the system can be specified if required; otherwise they are all assumed to be unity. 

		The magnitudes and coordinates can either be specified by the user or generated randomly. For now, coordinates must be specified in units of pixels; sub-pixel positioning of stars has not been implemented.
	"""
	detector_height_px, detector_width_px = detector_size_px[0:2]
	if N_stars != None:
		m = np.random.uniform(low=m_max, high=m_min, size=(N_stars))
		coords = np.zeros((2,N_stars))
		coords[0,:] = np.random.uniform(high=detector_height_px, size=(1, N_stars))
		coords[1,:] = np.random.uniform(high=detector_width_px, size=(1, N_stars))
	elif m == None or coords == None:
		print 'ERROR: ether both m and coords OR N_stars be specified!'
		return
	else:
		if not(np.isscalar(m)):
			N_stars = len(m)
			if len(m) != coords.shape[0]:
				print 'ERROR: the length of the star magnitudes vector must be equal to that of the first dimension of the coords array!'
				return
		else:
			N_stars = 1
			m_tmp = m 
			m = np.ndarray(shape=(1,1))
			m[0] = m_tmp
			coords_tmp = coords
			coords = np.ndarray(shape=(2,1))
			coords[:,0] = coords_tmp		

	if wavelength_m == None and band != None:
		wavelength_m = FILTER_BANDS_M[band][0]
		bandwidth_m = FILTER_BANDS_M[band][1]
	
	starfield = np.zeros(detector_size_px)
	for k in range(N_stars):
		# Total expected e/s from the star from the whole aperture (NOT per pixel)
		Sigma_electrons = surfaceBrightness2countRate(mu = m[k], A_tel = A_tel, tau = tau, qe = qe, gain = gain, magnitudeSystem = magnitudeSystem, wavelength_m = wavelength_m, bandwidth_m = bandwidth_m)

		# Adding the star's contribution to the image.
		star, I, P_0, P_sum, I_0 = airyDisc(wavelength_m = wavelength_m, f_ratio = f_ratio, l_px_m = l_px_m, detector_size_px = detector_size_px, coords = coords[:,k], P_0 = Sigma_electrons)

		# Multiplying by the exposure time to get the electron count.
		starfield += star * t_exp	

	if verbose:
		print '-------------------------------------------------------------------------------------------------------'
		print '#\tCoordinates\t\tMagnitude\t\tP_0\t\tP_sum\t% in image'
		print '-------------------------------------------------------------------------------------------------------'
		for k in range(N_stars):
			print '%d:\t(%6.2f, %6.2f)\t%4.2f\t\t%15.2f\t%15.2f\t%4.2f' % (k+1, coords[0,k], coords[1,k], m[k], P_0, P_sum, 100*P_sum/P_0)
		print '-------------------------------------------------------------------------------------------------------'

	# Converting to image counts
	image_count = expectedCount2count(starfield, detectorSaturation = detectorSaturation)

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE,FIGSIZE))
		if len(m) == 1:
			plt.suptitle(r'Starfield, $m = %.2f$' % m)
		else:
			plt.suptitle('Starfield')
		plt.subplot(1,2,1)
		plt.imshow(starfield, norm=LogNorm())
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Expected electron flux')
		plt.subplot(1,2,2)
		if max(starfield.flatten() > 0):
			plt.imshow(image_count, norm=LogNorm())
		else:
			plt.imshow(image_count)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Simulated image')
		plt.show()

	return image_count, starfield, m, coords