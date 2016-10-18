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

import miscutils as mu
import numpy as np
import pdb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary_r')
import os

# linguine modules 
from linguineglobals import *
import etcutils
import imutils

############################################################################################
# The following functions are for use with GALFIT.
############################################################################################
def simulateSersicGalaxy(im_out_fname, # Output FITS file name
	height_px, 	# Height of output file
	width_px, 	# Width of output file
	mu_e,		# Surface brightness magnitude at effective radius
	R_e_px,		# Effective (half-light) radius (can be float)
	n, 			# Sersic index
	plate_scale_as_px, 	# Plate scale
	axis_ratio,	# Axis ratio (a/b)
	zeropoint = -AB_MAGNITUDE_ZEROPOINT, # Careful of the minus sign!
	pos_px = None, # Position of galaxy in frame 
	PA_deg = 0,	# Rotation angle
	object_type = 'sersic2',
	galfit_input_fname = "galfit_input.txt",
	plotIt = False,
	overwriteExisting = False
	):
	"""
		Return a simulated image of a galaxy (made using GALFIT) given the inputs.
	"""
	
	if not im_out_fname.endswith('.fits'):
		im_out_fname += '.fits'

	# Writing the parameters file.
	if not os.path.isfile(im_out_fname) or overwriteExisting:		
		galfit_input_fname, im_out_fname = writeGALFITparamsFile(galfit_input_fname, im_out_fname, height_px, width_px, mu_e, R_e_px, n, plate_scale_as_px, axis_ratio, zeropoint, pos_px, PA_deg, object_type)
		# Calling GALFIT.
		callGALFIT(galfit_input_fname)
	else:
		print("WARNING: I found a GALFIT FITS file '{}' with the same name as the input filename, so I am using it instead!".format(galfit_input_fname))

	# Editing the header to include the input parameters.
	hdulist = astropy.io.fits.open(im_out_fname, mode='update')
	if overwriteExisting:
		hdulist[0].header['R_E_PX'] = R_e_px
		hdulist[0].header['MU_E'] = mu_e
		hdulist[0].header['SER_IDX'] = n
	im_raw = hdulist[0].data
	hdulist.flush()
	hdulist.close()

	# Plotting.
	if plotIt:
		mu.newfigure()
		plt.imshow(im_raw)
		plt.title("GALFIT-generated image")
		mu.colorbar()
		plt.show()

	return im_raw

############################################################################################
from subprocess import call
def callGALFIT(galfit_input_fname):
	""" 
		Call GALFIT on the file galfit_input_fname.
	"""
	print("Calling GALFIT...")
	call(["galfit", galfit_input_fname])

############################################################################################
def writeGALFITparamsFile(
	galfit_input_fname,	# Name of GALFIT input parameters file
	im_out_fname, 		# Output FITS file name
	height_px, 			# Height of output file
	width_px, 			# Width of output file
	mu_e,				# Surface brightness magnitude at effective radius
	R_e_px,				# Effective (half-light) radius (can be float)
	n, 					# Sersic index
	plate_scale_as_px, 	# Plate scale
	axis_ratio,			# Axis ratio (a/b)
	zeropoint = -AB_MAGNITUDE_ZEROPOINT, # Careful of the minus sign!
	pos_px = None, 		# Position of galaxy in frame 
	PA_deg = 0,	# Rotation angle
	object_type = 'sersic2',
	):
	"""
		Write a .txt file that can be used as input to GALFIT containing control and fitting parameters.
	"""
	
	# By default, the galaxy is centered in the middle of the image plane.
	if not pos_px:
		# pos_px = (height_px/2, width_px/2)
		pos_px = (width_px/2, height_px/2)

	# Checking the file names.
	if galfit_input_fname.lower().endswith('.txt') == False:
		galfit_input_fname += '.txt'
	if im_out_fname.lower().endswith('.fits') == False:
		im_out_fname += '.fits'

	# Writing the input parameters file.
	print("Writing GALFIT input parameters profile to file {}...".format(galfit_input_fname))
	with open(galfit_input_fname, "w") as text_file:
		text_file.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		text_file.write("B)\t{}\n".format(im_out_fname))
		# text_file.write("H)\t0\t{:d}\t0\t{:d}\n".format(height_px, width_px))
		text_file.write("H)\t0\t{:d}\t0\t{:d}\n".format(width_px, height_px))
		text_file.write("J)\t{:.10f}\n".format(zeropoint))
		text_file.write("K)\t{:.10f}\t{:.10f}\n".format(plate_scale_as_px, plate_scale_as_px))
		text_file.write("O)\tregular\n")
		text_file.write("P)\t0\n")
		# Object fitting
		text_file.write("\n# INITIAL FITTING PARAMETERS\n")
		text_file.write("# Object 1\n")
		text_file.write("0)\t{}\n".format(object_type))
		text_file.write("1)\t{:.10f}\t{:.10f}\n".format(pos_px[0], pos_px[1]))
		text_file.write("3)\t{:.10f}\n".format(mu_e))
		text_file.write("4)\t{:.10f}\n".format(R_e_px))
		text_file.write("5)\t{:d}\n".format(n))
		text_file.write("9)\t{:.10f}\n".format(axis_ratio))
		text_file.write("10)\t{:.10f}\n".format(PA_deg))
		text_file.write("Z)\t0\n")

	return galfit_input_fname, im_out_fname

############################################################################################
# The following functions are for generating an image of a galaxy without GALFIT (not recommended for now)
############################################################################################
def sersic(n, R_e, R, mu_e,
	zeropoint = 0,
	wavelength_m = None
	):
	"""
		Returns the Sersic profile at radial distance R in units of surface brightness mu
		(given in magnitudes/arcsec^2, by default in AB magnitudes) and in flux units 
		with Sersic index n given half-light radius R_e and mu(R=R_e) = mu_e.
	"""
	# Calculating the constraining parameter.
	F_e = etcutils.surface_brightness2flux(mu = mu_e, zeropoint=zeropoint, wavelength_m=wavelength_m)

	# Calculating b_n given the Sersic index n.
	if n > 0.5 and n < 8:
		b_n = 1.9992 * n - 0.3271
	elif n >= 8:
		b_n = 2 * n - 1 / 3
	else:
		print('ERROR: haven\'t implemented n < 0.5 yet! Check again later...')
		return

	F = {
		'F_nu_cgs' 		: F_e['F_nu_cgs'] * np.exp(- b_n * (np.power(R/R_e, 1/n) - 1)),
		'F_lambda_cgs' 	: F_e['F_lambda_cgs'] * np.exp(- b_n * (np.power(R/R_e, 1/n) - 1)) ,
		'F_nu_si' 		: F_e['F_nu_si'] * np.exp(- b_n * (np.power(R/R_e, 1/n) - 1)),
		'F_lambda_si' 	: F_e['F_lambda_si'] * np.exp(- b_n * (np.power(R/R_e, 1/n) - 1))
	}
	mu = - 2.5 * np.log10(F['F_nu_cgs']) - zeropoint

	return R, mu, F

############################################################################################
def sersic2D(n, R_e, mu_e,
	theta_rad = 0,		# Angle between major axis and detector horizontal (radians)
	i_rad = 0,			# Inclination angle (radians; face-on corresponds to i = 0)
	R_max = None,		# Plotting limit. Default is 20 * R_e
	R_trunc = np.inf, 	# Disc truncation radius
	gridsize = 500,		# Number of returned grid points
	zeropoint = 0,
	wavelength_m = None,
	plotIt = False,
	R_units = 'kpc'
	):	
	" Returns 2D Sersic intensity and surface brightness plots. "
	if R_max == None:
		R_max = 20 * R_e

	# Making a 2D intensity plot of the galaxy given its inclination and orientation.
	dR = 2 * R_max / gridsize
	scaleFactor = 2
	imsize = gridsize*scaleFactor
	r = np.linspace(-R_max*scaleFactor, +R_max*scaleFactor, imsize)
	X, Y = np.meshgrid(r, r)
	if theta_rad != 0:
		print('WARNING: rotations are still kinda dodgy. Proceed with caution!')
	R = imutils.rotateAndCrop(image_in_array = np.sqrt(X * X + Y * Y / (np.cos(i_rad) * np.cos(i_rad))), angle=theta_rad * 180 / np.pi, cropArg=(imsize-gridsize)//2)
	# Calculating the Sersic flux and surface brightness profiles
	R, mu_map, F_map = sersic(n=n, R_e=R_e, R=R, mu_e=mu_e, zeropoint=zeropoint, wavelength_m=wavelength_m)
	# Truncating the profiles
	mu_map[R>R_trunc] = np.inf
	for key in F_map:
		F_map[key][R>R_trunc] = 0

	if plotIt:
		mu.newfigure(2,1)
		plt.subplot(1,2,1)
		plt.imshow(F_map['F_nu_cgs'], norm=LogNorm(), extent = [-dR*gridsize/2,dR*gridsize/2,-dR*gridsize/2,dR*gridsize/2])
		plt.xlabel(r'$R$ (%s)' % R_units)
		plt.ylabel(r'$R$ (%s)' % R_units)
		mu.colorbar()
		plt.title('2D intensity map')
		plt.subplot(1,2,2)
		plt.imshow(mu_map, extent = [-dR*gridsize/2,dR*gridsize/2,-dR*gridsize/2,dR*gridsize/2])
		plt.xlabel(r'$R$ (%s)' % R_units)
		plt.ylabel(r'$R$ (%s)' % R_units)
		mu.colorbar()
		plt.title('2D surface brightness map')
		plt.suptitle('2D Sersic profiles')
		plt.show()

	return R, dR, F_map, mu_map

############################################################################################
def exportGalaxyFITSFile(image_in_array, n, R_e, mu_e, z, R_trunc, i_deg, band, seeing_as, t_exp, N_exp,
	overwriteExisting = True,
	relpath = None
	):
	""" 
		Export a FITS file to be used as input to GALFIT.
	"""
	headerData = {
		'EXPTIME' : t_exp,
		'NCOMBINE' : 1,
		'GAIN' : 1, 	# if the image counts are in units of electrons
		# 'GAIN' : detector.adu_gain,  # if the image counts are in units of ADU
		'RDNOISE' : np.sqrt(N_exp) * detector.read_noise * detector.read_noise
	}
	if R_trunc == np.inf:
		R_trunc = 0
	fname = 'gal_%d_%02d_%02d_%03d_%02d_%01d__%s_%01d_%03d_%03d' % (n, R_e, mu_e, z*1e3, R_trunc, i_deg, band, seeing_as, t_exp*1e3, N_exp)
	if relpath != None:
		fname = relpath + '/' + fname
	print("Exporting image data to file: ", fname,".fits")
	imutils.exportFitsFile(image_in_array = image_in_array, fname = fname, otherHeaderData = headerData, overwriteExisting = overwriteExisting)
