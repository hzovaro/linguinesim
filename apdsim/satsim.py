############################################################################################
#
# 	File:		galsim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	For simulating images of galaxies. 
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
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
############################################################################################
from __future__ import division, print_function
from linguinesim.apdsim import *

# In this module, we want to generate an image of a satellite.
# The image must be returned as im_raw to be consisted with galsim.py, to enable us to switch between simulating galaxies and satellites with ease.
# The image must be given in units of expected count rate.

"""
	~~* Brainthoughts *~~

	- Francis: satellite magnitudes ~4-5 (AB?): totally source-limited! So sky background calculations probably aren't too critical. Still want to estimate them (how? empirically?) and use them but they won't affect the result to the same extent that they have for the galaxy simulations.
	- Need to implement the stochastic gain process for eMCCDs in countrate2count(). Perhaps make a Detector class.
	- Maximum frame rate ~60 Hz so t_exp = 1/60 s minimum.
	- We need to look at LI simulation output as a function of 
		- altitude;
		- satellite geometry (shape and size);
		- reflectivity; 
		- time of day (sun position);
		- orbit (i.e. track in sky)
	- Francis: treat satellite as a Lambertian surface. 
	- TO START OFF WITH: 
		0. Super basic: treat as a point source with some magnitude. 
		1. Spherical object w/ Lambertian surface. Calculate magnitude as a function of phase angle for a sphere. 
		2. Images showing how magnitude changes with phase angle.
	- Write a function that calculates the phase angle given
		- Groundstation coordinates 
		- Object position in the sky - what is the best way to quantify this?
		- Time of day (i.e. sun position)

"""

def simulateSatellite(dim_m,
	height_px, 			# Height of output file
	width_px, 			# Width of output file
	alt_m,				# Altitude of satellite
	alpha_rad,			# Angle between surface normal and sun (radians)
	plate_scale_as_px,
	band,				# Imaging band
	A_collecting_m2, 	# Collecting area of telescope (m^2)
	t_exp,				# Exposure time (s)
	albedo=1.0,			# Albedo (fraction of reflected light)
	im_out_fname=None 	# Output FITS file name
	):
	"""

	Return an image of a satellite's solar panel (modelled as a Lambertian reflector) illuminated by the sun with a given phase angle.
	
	(A side note: some satellites will have sun-tracking solar panels!)
	(Another side note: simulating GEO satellites will be much easier as they will appear as an unresolved blob... so all we need is the magnitude and we're done.)

	Parameters
	-----------
	mag: float
		AB magnitude of the satellite.
	alt_m: float
		Altitude of the satellite in metres.
	alpha_rad: float
		Phase angle; that is, the angle between the surface normal and the line connecting the surface to the illuminating source (the sun).
	theta_rad: float
		Angle between the surface normal and the line connecting the surface to the observer. 
	shape: string
		Morphology of the satellite. For now, 'box' and 'sphere' are valid entries.
	dim_m: tuple
		Dimensions of the satellite in metres. 
		shape == 'box': dim_m = (height_m, width_m, theta_deg)
		shape == 'sphere': dim_m = (radius)
	plate_scale_as_px: float
		Plate scale of the detector in arcseconds per pixel.
	height_px: int
		Height of output image.
	width_px: int
		Width of output image.
	(optional) im_out_fname: string
		If specified, will save the image to a FITS file with filename im_out_fname.
	"""


	# 'Truth' image will be undersampled. So need to convolve it with a PSF prior to rotating it. As such the PSF needs to be counter-rotated by the same amount (not so important for now).

	# To begin with, let's simulate a solar panel. So basically a big, flat rectangle. 
	# 0. Convert the physical dimensions of the object to arcseconds so we can make an image of it on our detector.
	dim_as = tuple(3600 * np.rad2deg(l_m / alt_m) for l_m in dim_m)
	dim_px = tuple(np.round(l_as / plate_scale_as_px) for l_as in dim_as)
	A_as = dim_as[0] * dim_as[1]

	# Make the image.
	x = np.arange(height_px)-height_px/2
	y = np.arange(width_px)-width_px/2
	X, Y = np.meshgrid(y,x)
	idxs = np.where( (Y < dim_px[0]/2) * (Y > (-dim_px[0]/2)) * (X < dim_px[1]/2) * (X > (-dim_px[1]/2)))
	im_raw = np.zeros( (height_px, width_px) )

	# Determining the surface brightness and magnitude. 
	I_si = solarFlux(band) * albedo * np.cos(alpha_rad) / 2 / np.pi 			# Irradiance of the emitting surface (W/m^2/sr). Only depends on angle between surface normal and illuminating source (the sun)
	I_cgs = I_si * 1e7 * 1e-4 / ((3600 * 180 / np.pi)**2)						# ergs/s/cm^2/arcsec^2
	F_lambda_cgs = I_cgs / FILTER_BANDS_M[band][1]								# ergs/s/cm^2/m/arcsec^2
	F_nu_cgs = F_lambda_cgs * (FILTER_BANDS_M[band][0]**2) / constants.c 		# ergs/s/cm^2/Hz/arcsec^2
	mu = - 2.5 * np.log10(F_nu_cgs) - AB_MAGNITUDE_ZEROPOINT 					# AB mag per arcsec^2
	mag = - 2.5 * np.log10(F_nu_cgs * A_as) - AB_MAGNITUDE_ZEROPOINT			# AB mag
	# So regardless of the viewing angle theta of the surface, the surface brightness is constant. But the TOTAL flux we receive from the source depends on A*cos(viewing angle), so the INTEGRATED magnitude (surface brightess * arcsec^2 of solar panel) will also change.

	# Scaling to a magnitude: multiply up by the photon count.
	E_phot = constants.h * constants.c / FILTER_BANDS_M[band][0]	# Photon energy (W)
	im_raw[idxs] = I_si * (np.deg2rad(plate_scale_as_px / 3600))**2 * A_collecting_m2 / E_phot * t_exp	# photons/pixel
	# im_raw *= tau * qe * gain	# counts/pixel

	# Saving to file
	if im_out_fname:
		exportFitsFile(image_in_array=im_raw, fname=fname, otherHeaderData={
			'dim_m' : dim_m,
			'height_px' : height_px, 			# Height of output file
			'width_px': width_px, 			# Width of output file
			'alt_m':alt_m,				# Altitude of satellite
			'alpha_rad':alpha_rad,			# Angle between surface normal and sun (radians)
			'plate_scale_as_px':plate_scale_as_px,
			'band':band,				# Imaging band
			'A_collecting_m2':A_collecting_m2, 	# Collecting area of telescope (m^2)
			't_exp':t_exp				# Exposure time (s)
			})

	return im_raw, mu, mag

def solarFlux(band):
	"""
	Calculate the total solar flux (in W/m^2) integrated over the given band by approximating the sun to a black body with a temperature of 5777 K.

	Parameters
	-----------
	band: dict key
		Band over which to compute the integrated solar flux.

	Returns:
	-----------
		The integrated solar flux in the specified band in units of W/m^2.
	"""
	E_phot = constants.h * constants. c / FILTER_BANDS_M[band][0]	# We must multiply the output by E_phot because thermalEmissionIntensity returns the irradiance in photons/s, not in W
	return thermalEmissionIntensity(T=T_SUN_K,
		wavelength_min=FILTER_BANDS_M[band][2],
		wavelength_max=FILTER_BANDS_M[band][3],
		Omega=(2 * R_SUN_M / DIST_SUN_M)**2, # Solid angle subtended by sun at the observer's position
		A=1) * E_phot		# Collecting area (so the cross-sectional area of the satellite). Because we assume a Lambertian surface, we define A here to be the area of the solar panel, and NOT the area projected along the line-of-sight from the solar panel to the sun (that comes lateR). Set as 1 m^2 so we get the answer back in units of W/m^2.