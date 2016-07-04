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
#	- Better formula for calculating b_n.
#	- Units of I_e: assume that mu_e is expressed in units of AB magnitudes/arcsec^2.
#
############################################################################################
from __future__ import division
from apdsim import *

def sersic(n, R_e, R, mu_e,
	zeropoint = 0,
	wavelength_m = None
	):
	"""
		Returns the Sersic profile at radial distance R in units of surface brightness mu
		(given in magnitudes/arcsec^2, by fault in AB magnitudes) and in flux units 
		with Sersic index n given half-light radius R_e and mu(R=R_e) = mu_e.
	"""
	# Calculating the constraining parameter.
	F_e = surfaceBrightnessToFlux(mu = mu_e, zeropoint=zeropoint, wavelength_m=wavelength_m)

	# Calculating b_n given the Sersic index n.
	if n > 0.5 and n < 8:
		b_n = 1.9992 * n - 0.3271
	elif n >= 8:
		b_n = 2 * n - 1 / 3
	else:
		print 'ERROR: haven\'t implemented n < 0.5 yet! Check again later...'
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
		print 'WARNING: rotations are still kinda dodgy. Proceed with caution!'
	R = rotateAndCrop(image_in_array = np.sqrt(X * X + Y * Y / (np.cos(i_rad) * np.cos(i_rad))), angle=theta_rad * 180 / np.pi, cropArg=(imsize-gridsize)//2)
	# Calculating the Sersic flux and surface brightness profiles
	R, mu_map, F_map = sersic(n=n, R_e=R_e, R=R, mu_e=mu_e, zeropoint=zeropoint, wavelength_m=wavelength_m)
	# Truncating the profiles
	mu_map[R>R_trunc] = np.inf
	for key in F_map:
		F_map[key][R>R_trunc] = 0

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE,FIGSIZE))
		plt.subplot(1,2,1)
		plt.imshow(F_map['F_nu_cgs'], norm=LogNorm(), extent = [-dR*gridsize/2,dR*gridsize/2,-dR*gridsize/2,dR*gridsize/2])
		plt.xlabel(r'$R$ (%s)' % R_units)
		plt.ylabel(r'$R$ (%s)' % R_units)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('2D intensity map')
		plt.subplot(1,2,2)
		plt.imshow(mu_map, extent = [-dR*gridsize/2,dR*gridsize/2,-dR*gridsize/2,dR*gridsize/2])
		plt.xlabel(r'$R$ (%s)' % R_units)
		plt.ylabel(r'$R$ (%s)' % R_units)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('2D surface brightness map')
		plt.suptitle('2D Sersic profiles')
		plt.show()

	return R, dR, F_map, mu_map

############################################################################################
def exportGalaxyFITSFile(image_in_array, n, R_e, mu_e, z, R_trunc, i_deg, band, seeing_as, t_exp, N_exp,
	overwriteExisting = True):
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
	fname = 'gal_%d_%02d_%02d_%03d_%02d_%01d__%s_%01d_%03d_%03d' % (n, R_e, mu_e, z*1e3, R_trunc, i_deg, band, seeing_as, t_exp*1e3, N_exp)
	print "Exporting image data to file: ", fname, ".fits"
	exportFITSFile(image_in_array = image_in_array, fname = fname, otherHeaderData = headerData, overwriteExisting = overwriteExisting)


