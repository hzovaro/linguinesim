############################################################################################
#
# 	File:		imutils.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		29/06/2016
#
#	Description:
#	Image processing utilities.
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
from __future__ import division
from apdsim import *

"""
	Handy notes for FITS files:
		- to see the header information: print(repr(hdulist[0].header))
		- exposure time: hdulist[0].header['EXPTIME']
		- Other header keys:
			- ATODGAIN, BANDWID, CENTRWV, INSTRUME, TELESCOP, OBJECT, NAXIS, NAXIS1, NAXIS2,
"""
def getRawImages(fname, 
	plotIt=False,
	idx=0):
	" Return an array of the image(s) stored in the FITS file fname. "
	hdulist = fits.open(fname)
	images_raw = hdulist[idx].data
	hdulist.close()

	images_raw, N, height, width = getImageSize(images_raw)
	if plotIt:
		for k in range(N):
			plt.figure()
			plt.imshow(images_raw[k])
			plt.title('Raw image %d from FITS file' % k + 1)
		plt.show()

	return np.squeeze(images_raw), hdulist

#########################################################################################################
def imageToArray(im):
	" Convert an Image object im into an array. "
	width, height = im.size
	image_map = list(im.getdata())
	image_map = np.array(image_map).reshape(height, width)
	return image_map

#########################################################################################################
def rotateAndCrop(image_in_array, angle, cropArg, 
	plotIt=False):
	" Rotate and crop an array of N images stored in ndarray image_in_array counterclockwise by a given angle and then crop the image using coordinates (left, upper, right, lower) "
	image_in_array, N, height, width = getImageSize(image_in_array)

	# Crop options:
	#	1. One number given: crop by the same amount on all sides
	if type(cropArg) == int:
		cropIdx = (cropArg, cropArg, width - cropArg, height - cropArg)
	#	2. Two numbers given: crop by the same width and height either side.
	elif type(cropArg) == tuple and len(cropArg) == 2:
		cropIdx = (cropArg[1], cropArg[0], width - cropArg[1], height - cropArg[0])
	#	3. Four numbers given: crop input is given as the (left, top, right, bottom) indices.
	elif type(cropArg) == tuple and len(cropArg) == 4:
		cropIdx = cropArg

	# Convert to an Image object.
	for k in range(N):
		image = Image.fromarray(image_in_array[k])
		# Rotate.
		if angle != 0.:
			image = image.rotate(angle)
		# Crop.
		image = image.crop(cropIdx)
		image = imageToArray(image)
		if k == 0:
			image_out_array = np.ndarray((N, image.shape[0], image.shape[1]))
		# Convert back to an array.
		image_out_array[k] = image

	if plotIt:
		plt.figure()
		plt.subplot(1,2,1)
		plt.imshow(image_in_array[0])
		plt.title('Input image')
		plt.subplot(1,2,2)
		plt.imshow(image_out_array[0])
		plt.title('Output image')
		plt.suptitle('Rotating and cropping image')
		plt.show()

	return np.squeeze(image_out_array)

#########################################################################################################
def getImageSize(image_in_array):
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
		print "ERROR: invalid image array shape!"
		return -1
	return (image_out_array, N, height, width)