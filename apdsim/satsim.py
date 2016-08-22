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
from __future__ import division
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
		1. Spherical object w/ Lambertian surface.
		2. Images showing how magnitude changes with phase angle.
	- Write a function that calculates the phase angle given
		- Groundstation coordinates 
		- Object position in the sky - what is the best way to quantify this?
		- Time of day (i.e. sun position)

"""

def simulateSatellite(mag, alt_m, shape, dim_m, plate_scale_as_px,
	height_px, 			# Height of output file
	width_px, 			# Width of output file
	im_out_fname=None 	# Output FITS file name
	):
	"""
	Parameters
	-----------
	mag: float
		AB magnitude of the satellite.
	alt_m: float
		Altitude of the satellite in metres.
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
	im_raw = np.zeros( (height_px, width_px) )

	# Scaling to a magnitude: normalise the image so that the sum of all pixel values are 1. Then, multiply up by the magnitude.
	im_raw /= 

	return im_raw