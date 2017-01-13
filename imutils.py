################################################################################
#
# 	File:		imutils.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Image processing utilities.
#
#	Copyright (C) 2016 Anna Zovaro
#
################################################################################
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
################################################################################
from __future__ import division, print_function
import miscutils as mu
import numpy as np
import ipdb
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary_r')
import pyfftw
import fftwconvolve

import astropy.io.fits

# linguine modules 
from linguineglobals import *

"""
	Handy notes for FITS files:
		- to see the header information: print(repr(hdulist[0].header))
		- exposure time: hdulist[0].header['EXPTIME']
		- Other header keys:
			- ATODGAIN, BANDWID, CENTRWV, INSTRUME, TELESCOP, OBJECT, NAXIS, NAXIS1, NAXIS2,
"""
def image_from_fits(fname, 
	plotit=False,
	idx=0):
	" Return an array of the image(s) stored in the FITS file fname. "
	if not fname.lower().endswith('fits'):
		fname += '.fits'
	hdulist = astropy.io.fits.open(fname)
	images_raw = hdulist[idx].data
	hdulist.close()

	images_raw, N, height, width = get_image_size(images_raw)
	if plotit:
		for k in range(N):
			plt.imshow(images_raw[k])
			plt.title('Raw image %d from FITS file' % (k + 1))
			plt.pause(0.1)
		plt.show()

	return np.squeeze(images_raw), hdulist

################################################################################
def image_obj_to_array(im):
	" Convert an Image object im into an array. "
	width, height = im.size
	image_map = list(im.getdata())
	image_map = np.array(image_map).reshape(height, width)
	return image_map

################################################################################
def centre_crop(im, sz_final, 
		units = 'px',
		plate_scale_as_px = 1,
		centre_coords_rel = np.array([0,0])):
	"""
		Crop an array by equal amounts on each side. 

		If desired, the input units can be specified in arcsec. 
	"""

	if units=='arcsec':
		if type(sz_final) == tuple:
			sz_final = np.array(sz_final)
		if type(centre_coords_rel) == tuple:
			centre_coords_rel = np.array(centre_coords_rel)
		sz_final[0] = int(np.round(sz_final[0] / plate_scale_as_px))
		sz_final[1] = int(np.round(sz_final[1] / plate_scale_as_px))
		centre_coords_rel[0] = int(np.round(centre_coords_rel[0] / plate_scale_as_px))
		centre_coords_rel[1] = int(np.round(centre_coords_rel[1] / plate_scale_as_px))

	im = get_image_size(im)[0]

	if matplotlib.cbook.is_scalar(sz_final):
		sz_height = sz_final
		sz_width = sz_final
	else:
		sz_height = sz_final[0]
		sz_width = sz_final[1]

	crop_height = max((im.shape[1] - sz_height) // 2, 0)
	crop_width = max((im.shape[2] - sz_width) // 2, 0)
	
	im_cropped = im[
		:,
		crop_height + centre_coords_rel[0] : crop_height + min(im.shape[1], sz_height) + centre_coords_rel[0], 
		crop_width  + centre_coords_rel[1] : crop_width + min(im.shape[2], sz_width) + centre_coords_rel[1] 
		]

	return np.squeeze(im_cropped)

################################################################################
def get_image_size(image_in_array):
	" This function takes as input an image array which is either 2- or 3-dimensional. It returns an new array with dimensions (N, height, width). This function is basically to ensure consistency in this module in how images are stored (to eliminate the need to constantly check the dimensions of input image arguments) "
	if len(image_in_array.shape) == 3:
		N = image_in_array.shape[0]
		height = image_in_array.shape[1]
		width = image_in_array.shape[2]
		image_out_array = image_in_array	
	elif len(image_in_array.shape) == 2:
		N = 1
		height = image_in_array.shape[0]
		width = image_in_array.shape[1]
		image_out_array = np.ndarray((1, height, width), dtype=type(image_in_array[0,0]))
		image_out_array[0,:,:] = image_in_array
	else:
		print("ERROR: invalid image array shape!")
		return -1
	return (image_out_array, N, height, width)

################################################################################
def export_fits(image_in_array, fname,
	otherHeaderData = None, 
	overwrite_existing = False):
	"""
		Exports a FITS file containing the image data contained in image_in_array and optional header 
		data stored in the dictionary otherHeaderData. Does not overwrite existing file fname.fits 
		by default, but will if overwrite_existing is set to True.

		Some notes: 
		- It is generally better to give the pixel values in terms of ADU count units instead of flux 
		values as the sky background component usually appears deceptively small
		- If the image is an averaged image of, say, N exposures, then NCOMBINE = N and GAIN and
		RDNOISE are those values corresponding to a single exposure. 
		- RDNOISE is expressed in units of electrons (not electrons rms?)
		- If the image is a stacked (summed) image, then NCOMBINE = 1, GAIN corresponds to that of 
		a single image but RDNOISE = sqrt(N) * RN where RN corresponds to that of a single image.

		Note: HDU stands for 'Header Data Unit'
	"""
	hdu = astropy.io.fits.PrimaryHDU()

	# Add data. Note that this automatically updates the header data with the axis sizes. 
	hdu.data = image_in_array

	# Optional header data:
	if otherHeaderData != None:
		for key in otherHeaderData:
			hdu.header[key] = otherHeaderData[key]
		
	# Write to disk.
	if fname.endswith('.fits') == False:
		fname += '.fits'
	hdu.writeto(fname, clobber = overwrite_existing)


################################################################################
def fourier_resize(im, scale_factor,
	conserve_pixel_sum=True):

	# Resize an image using a Fourier transform.
	h,w = im.shape
	h_s = h/scale_factor
	w_s = w/scale_factor

	# Take the Fourier transform & shift so that low frequency components are in
	# the centre.
	im_fft = np.fft.fftshift(pyfftw.interfaces.numpy_fft.fft2(im))

	# Crop it.
	im_fft_cropped = centre_crop(im_fft,sz_final=(h_s,w_s))

	# Inverse transform.
	im_resized = np.abs(pyfftw.interfaces.numpy_fft.ifft2(
		np.fft.fftshift(im_fft_cropped)))

	if conserve_pixel_sum:
		sum_before = sum(im.flatten())
		sum_after = sum(im_resized.flatten())
		im_resized = im_resized / sum_after * sum_before

	return im_resized

################################################################################
def gaussian_smooth(im, sigma):
	# Smooth an image by convolving it with a Gaussian kernel.
	# Generate a Gaussian kernel.
	# Make it extend out to 5sigma.
	x = np.arange(-5*sigma, +5*sigma, dtype='float')
	y = x
	X, Y = np.meshgrid(x,y)
	kernel = 1 / (2 * np.pi * sigma**2) * \
	np.exp( -(X**2 + Y**2) / (2 * sigma**2) )
	return fftwconvolve.fftconvolve(im, kernel, mode='same')
