##########################################################################################################
#
# 	File:		lisim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A module for simulating lucky imaging.
#
#	Copyright (C) 2016 Anna Zovaro
#
#	Lucky imaging techniques to implement:
#	- Shifting-and-stacking via
#		- Cross-correlation
#		- Aligning to brightest pixel
#		- Drizzle algorithm for image alignment (sub-integer alignment)
#
#	- Frame selection techniques:
#		- Rank in order of brightest pixel value
#		- Cross-correlate the ideal PSF (say, Airy disc) with a subsection of the image containing a guide star--peak of the x-corr indicates the correlation (basically the Strehl) whilst its position gives the shift that needs to be applied 
#		- Rank in order of the fraction of light concentrated in the brightest pixel of the guide star PSF
#
#########################################################################################################
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
####################################################################################################
from __future__ import division
from __future__ import print_function
from apdsim import *

def addTurbulence(images, N_tt, sigma_tt_px,
	crop_tt=None):
	""" 
		Add turbulence to an input `truth' image. Returns N_tt copies of the input image with randomised turbulence added. 
		Just tip and tilt for now with a standard deviation of sigma_tt_px in both dimensions.
	"""
	print("Adding randomised tip/tilt to image(s)", end="")

	# Tip and tilt for now	
	images, N, height, width = getImageSize(images)

	# Output array of images
	if crop_tt == None:
		images_tt = np.ndarray((N, N_tt, height, width))
	else:
		# if type(crop_tt) == int:
		if len(crop_tt.shape) == 0:
			images_tt = np.ndarray((N, N_tt, height - 2 * crop_tt, width - 2 * crop_tt))	
		else:
			images_tt = np.ndarray((N, N_tt, height - 2 * crop_tt[0], width - 2 * crop_tt[1]))

	# Array to hold the tip/tilt offsets
	tt_idxs = np.ndarray((N, N_tt, 2))
	
	for k in range(N):
		for j in range(N_tt):
			print('.', end="")
			shift_height = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
			shift_width = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
			image_tt = shift(images[k], (shift_height, shift_width))
			tt_idxs[k,j] = [shift_height, shift_width]

			# Cropping the image if necessary
			if crop_tt == None:
				images_tt[k,j] = image_tt
			else:
				images_tt[k,j] = rotateAndCrop(image_tt, angle=0., cropArg=crop_tt)
		print('\n')
	return np.squeeze(images_tt), np.squeeze(tt_idxs)

####################################################################################################
def shiftAndStack(images, 
	image_ref=None,
	N=None, 
	showAnimatedPlots=False,	# Show an animated plot window of the shifting-and-stacking process.
	plotIt=False):
	""" 
		Shift and stack the images given in the 3-dimensional array of N images. 
	"""
	print("Shifting and stacking images", end="")

	# Error checking
	# Need to convert to float if necessary.
	if type(images.flatten()[0]) != np.float64:
		images = images.astype(np.float64)
	if image_ref != None and type(image_ref.flatten()[0]) != np.float64:
		image_ref = image_ref.astype(np.float64)

	if len(images.shape) > 4:
		print("WARNING: for now, please only input a 3D array of images to shift and stack! I'm operating on the first set of images...")
		images = np.squeeze(images[0])	
	
	if len(images.shape) == 3:
		if N > images.shape[0]:
			print("ERROR: if specified, N must be equal to or less than the length of the first dimension of the images array.")
			return -1
		if image_ref == None:
			# If the reference image is not specified, we use the first image in the array as the reference: 
			# i.e. we align all other images to images[0].
			if N == None:
				N = images.shape[0]-1
			images = np.copy(images[1:])	# Only need to go through images 1:N-1.
			image_ref = np.copy(images[0])
		else:
			if N == None:
				N = images.shape[0]			
	else:
		# Error: cannot shift and stack a single image!
		print("ERROR: cannot shift and stack a single image! Input array must have N > 1.")
		return -1
	
	image_stacked = np.copy(image_ref)		# shifted-and-stacked image

	height = images[0].shape[0]
	width = images[0].shape[1]
	
	corrs = np.ndarray((N, 2 * height - 1, 2 * width - 1))	# Array to hold x-correlation results
	corr_peak_idxs = np.ndarray((N, 2))		# indices of the peak value in the x-correlation
	img_peak_idxs = np.ndarray((N, 2))		# shift in x and y computed from the x-correlation

	for k in range(N):
		print('.', end="")
		# Cross-correlate image k with the reference image to find the tip and tilt.
		corrs[k] = signal.fftconvolve(image_ref, images[k][::-1,::-1])
		corr_peak_idxs[k] = np.unravel_index(np.argmax(corrs[k]), (2 * height - 1, 2 * width - 1))
		img_peak_idxs[k][0] = - corr_peak_idxs[k][0] + (height - 1)
		img_peak_idxs[k][1] = - corr_peak_idxs[k][1] + (width - 1)

		# Shift-and-stack the images.
		# image_stacked += np.roll(np.roll(images[k], -img_peak_idxs[k][0].astype(int), 0), -img_peak_idxs[k][1].astype(int), 1)	
		image_stacked += shift(images[k], (-img_peak_idxs[k][0].astype(int), -img_peak_idxs[k][1].astype(int)))

		# Plotting
		if showAnimatedPlots:
			if k == 0:
				plt.figure(figsize=(3*FIGSIZE, FIGSIZE))
				plt.subplot(1,3,1)
				plt.imshow(images[0],origin='lower')
				plt.subplot(1,3,2)
				scat2 = plt.scatter(0.0,0.0,c='r',s=20)
				plt.subplot(1,3,3)
				scat3 = plt.scatter(0.0,0.0,c='g',s=40)

			plt.subplot(1,3,1)
			plt.imshow(image_stacked,origin='lower')

			plt.subplot(1,3,2)
			plt.imshow(images[k],origin='lower')	
			plotcoords = np.ndarray((2))
			plotcoords[1] = img_peak_idxs[k,0] + width / 2
			plotcoords[0] = img_peak_idxs[k,1] + height / 2
			scat2.set_offsets(plotcoords)

			plt.subplot(1,3,3)
			plt.imshow(corrs[k],interpolation='nearest',origin='lower')
			corr_peak_coords = np.ndarray((2))
			corr_peak_coords[0] = corr_peak_idxs[k][1]
			corr_peak_coords[1] = corr_peak_idxs[k][0]
			scat3.set_offsets(corr_peak_coords)

			# plt.scatter([peak_idxs[k][0]], [peak_idxs[k][1]], c='r', s=20)
			plt.draw()
			plt.pause(1)
	print('\n')

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE,FIGSIZE))
		plt.suptitle('Lucky imaging technique')
		plt.subplot(1,2,1)
		# plt.imshow(images[0], vmin=min(images[0].flatten()), vmax=max(image_stacked.flatten()))
		plt.imshow(images[0], norm=LogNorm())
		plt.title('Single exposure')
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.subplot(1,2,2)
		# plt.imshow(image_stacked, vmin=min(images[0].flatten()), vmax=max(image_stacked.flatten()))
		plt.imshow(image_stacked, norm=LogNorm())
		plt.title('Shifted-and-stacked image')
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.show()

	return image_stacked, img_peak_idxs

####################################################################################################
def luckyImaging(images, type, vararg):
	""" 
		Apply a Lucky Imaging (LI) technique to a sequence of images stored in the input array images. 
		The type of LI technique used is specified by input string type and any additional arguments which may be required are given in vararg.
	"""
	# Converting to an image array
	images, N, height, width = getImageSize(images)



	return image_lucky

####################################################################################################
def alignmentError(in_idxs, out_idxs,
	verbose=True):
	"""
		Compute the alignment errors arising in the Lucky Imaging shifting-and-stacking process given an input array of tip and tilt coordinates applied to the input images and the coordinates of the shifts applied in the shifting-and-stacking process. 
	"""
	N = in_idxs.shape[0]
	errs = np.zeros((N))
	n_errs = 0
		
	for k in range(N):
		errs[k] = np.sqrt(np.power(in_idxs[k,0] - out_idxs[k,0],2) + np.power(in_idxs[k,1] - out_idxs[k,1],2))
		if errs[k] > 0:
			n_errs += 1
			
	if verbose:
		print('------------------------------------------------')
		print('Tip/tilt coordinates\nInput\t\tOutput\t\tError')
		print('------------------------------------------------')
		for k in range(N):
			print('(%6.2f,%6.2f)\t(%6.2f,%6.2f)\t%4.2f' % (in_idxs[k,0],in_idxs[k,1],out_idxs[k,0],out_idxs[k,1],errs[k]))
		print('------------------------------------------------')

	return n_errs, errs
