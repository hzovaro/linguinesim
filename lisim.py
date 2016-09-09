####################################################################################################
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
####################################################################################################
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
####################################################################################################
from __future__ import division, print_function
import miscutils as mu
import numpy as np
import pdb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary_r')

import scipy.signal
import scipy.ndimage.interpolation
import astropy.modeling

# Multithreading/processing packages
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool	# dummy = Threads
from multiprocessing import Pool as ProcPool			# no dummy = Processes
import time

# linguine modules 
from linguineglobals import *
import fftwconvolve, obssim, etcutils, imutils


####################################################################################################
def luckyImage(im, psf, tt, noise_frame, scale_factor, etc_input):
	""" 
		This function can be used to generate a short-exposure 'lucky' image that can be input to the Lucky Imaging algorithms.
			Input: 	one 'raw' countrate image of a galaxy; one PSF with which to convolve it (at the same plate scale)
			Output: a 'Lucky' exposure. 			
			Process: convolve with PSF --> resize to detector --> add tip and tilt (from a premade vector of tip/tilt values) --> convert to counts --> add noise --> subtract the master sky/dark current. 
	"""
	# Convolve with PSF.
	im = obssim.convolvePSF(im, psf)		
	# Resize to detector.
	im = obssim.resizeImageToDetector(image_raw = im, source_plate_scale_as = 1, dest_plate_scale_as = scale_factor)
	# Add tip and tilt.
	im = obssim.addTipTilt_single(image = im, tt_idxs = tt)[0]
	# Convert to counts.
	im = etcutils.expectedCount2count(im, t_exp = etc_input['t_exp'])
	# Add noise. 
	im += noise_frame
	# Subtract the master sky/dark current. 
	# im = im.astype(np.float64) - master_frame
	return im

####################################################################################################
def shift_pp(image, img_ref_peak_idx, fsr, bid_area):
	if type(image) == list:
		image = np.array(image)	

	# Search in the bid area of the input image for the peak pixel coordinates.
	if bid_area:
		sub_image = imutils.rotateAndCrop(image_in_array = images, cropArg = bid_area)
	else:
		sub_image = image		
	img_peak_idx = np.asarray(np.unravel_index(np.argmax(sub_image), sub_image.shape))	

	# Shift the image by the relative amount.
	rel_shift_idx = (img_ref_peak_idx - img_peak_idx)
	image_shifted = scipy.ndimage.interpolation.shift(image, rel_shift_idx)

	peak_pixel_val = max(sub_image.flatten())	# Maximum pixel value (for now, not used)

	return image_shifted, -rel_shift_idx, peak_pixel_val

####################################################################################################
def shift_centroid(image, img_ref_peak_idx):
	if type(image) == list:
		image = np.array(image)

	img_peak_idx = _centroid(image)

	# Shift the image by the relative amount.
	rel_shift_idx = (img_ref_peak_idx - img_peak_idx)
	image_shifted = scipy.ndimage.interpolation.shift(image, rel_shift_idx)

	return image_shifted, -rel_shift_idx

####################################################################################################
def shift_xcorr(image, image_ref, buff, subPixelShift):
	if type(image) == list:
		image = np.array(image)

	height, width = image.shape
	if fftwconvolve.NTHREADS==0:
		corr = scipy.signal.fftconvolve(image_ref, image[::-1,::-1], 'same')
	else:
		corr = fftwconvolve.fftconvolve(image_ref, image[::-1,::-1], 'same')
	corr /= max(corr.flatten())	# The fitting here does not work if the pixels have large values!
	
	if subPixelShift: 
		# Fitting a Gaussian.
		Y, X = np.mgrid[-(height-2*buff)/2:(height-2*buff)/2, -(width-2*buff)/2:(width-2*buff)/2]
		try:		
			p_init = astropy.modeling.models.Gaussian2D(x_stddev=1.,y_stddev=1.)
		except:
			p_init = astropy.modeling.models.Gaussian2D(x_mean=1.,y_mean=1.,x_stddev=1.,y_stddev=1.,amplitue=1.)
		fit_p = astropy.modeling.fitting.LevMarLSQFitter()
		p_fit = fit_p(p_init, X, Y, corr[buff:height-buff, buff:width-buff])		
		rel_shift_idx = (p_fit.y_mean.value, p_fit.x_mean.value)	# NOTE: the indices have to be swapped around here for some reason!		
	else:
		rel_shift_idx = np.unravel_index(np.argmax(corr), corr.shape)
		rel_shift_idx = (rel_shift_idx[0] - height/2, rel_shift_idx[1] - width/2)
	
	image_shifted = scipy.ndimage.interpolation.shift(image, rel_shift_idx)	

	return image_shifted, tuple(-x for x in rel_shift_idx)

####################################################################################################
def luckyImaging(images, li_method, mode,
	image_ref = None,	# reference image
	fsr = 1,			# for peak pixel method
	bid_area = None,	# for peak pixel method
	N = None,
	subPixelShift = True,	# for xcorr method
	buff = 25, 			# for xcorr method
	timeIt = True
	):
	""" 
		Apply a Lucky Imaging (LI) technique to a sequence of images stored in the input array images. 
		The type of LI technique used is specified by input string type and any additional arguments which may be required are given in vararg.
	"""
	images, image_ref, N = _li_error_check(images, image_ref, N)
	if not timeIt:
		print("Applying Lucky Imaging technique '{}' to input series of {:d} images...".format(li_method, N))
	
	# For each of these functions, the output must be of the form 
	#	image_shifted, rel_shift_idxs	
	if li_method == 'xcorr':
		shift_fun = partial(shift_xcorr, image_ref=image_ref, buff=buff, subPixelShift=subPixelShift)
	
	elif li_method == 'peak_pixel':
		# Determining the reference coordinates.
		if bid_area:			
			sub_image_ref = imutils.rotateAndCrop(image_in_array = image_ref, cropArg = bid_area)
		else:
			sub_image_ref = image_ref
		img_ref_peak_idx = np.asarray(np.unravel_index(np.argmax(sub_image_ref), sub_image_ref.shape)) 

		shift_fun = partial(shift_pp, img_ref_peak_idx=img_ref_peak_idx, bid_area=bid_area, fsr=fsr)
	
	elif li_method == 'centroid':
		img_ref_peak_idx = _centroid(image_ref)
		shift_fun = partial(shift_centroid, img_ref_peak_idx=img_ref_peak_idx)
	
	else:
		print("ERROR: invalid Lucky Imaging method specified; must be 'xcorr', 'peak_pixel' or 'centroid' for now...")
		raise UserWarning

	# In here, want to parallelise the processing for *each image*. So make shift functions that work on a single image and return the shifted image, then stack it out here.
	
	tic = time.time()
	if mode == 'parallel':
		# Setting up to execute in parallel.
		images = images.tolist()	# Need to convert the image array to a list.

		# Executing in parallel.
		if NTHREADS == 0:
			# We are safe to use threads if NTHREADS==0 since we don't use pyfftw in this case.
			pool = ThreadPool()
		else:
			# Otherwise, we use processes instead.
			pool = ProcPool()
		results = pool.map(shift_fun, images, 1)
		pool.close()
		pool.join()

		# Extracting the output arguments.
		images_shifted = np.array(zip(*results)[0]) 
		rel_shift_idxs = np.array(zip(*results)[1])
		if li_method == 'peak_pixel' and fsr < 1:
			peak_pixel_vals = np.array(zip(*results)[2])

	elif mode == 'serial':
		# Loop through each image individually.
		images_shifted = np.zeros( (N, image_ref.shape[0], image_ref.shape[1]) )	
		rel_shift_idxs = np.zeros( (N, 2) )
		for k in range(N):
			if li_method == 'peak_pixel':
				if k == 0:
					peak_pixel_vals = np.zeros(N)
				images_shifted[k], rel_shift_idxs[k], peak_pixel_vals[k] = shift_fun(image=images[k])
			else:
				images_shifted[k], rel_shift_idxs[k] = shift_fun(image=images[k])
	else:
		print("ERROR: mode must be either parallel or serial!")
		raise UserWarning

	toc = time.time()
	if timeIt:
		print("APPLYING LUCKY IMAGING TECHNIQUE {}: Elapsed time for {:d} {}-by-{} images in {} mode: {:.5f}".format(li_method, N, image_ref.shape[0], image_ref.shape[1], mode, (toc-tic)))

	# If we're using an FSR < 1 in the peak pixel method, then we must do the following:
	#	1. Get our method to return a list of peak pixel values.
	#	2. Sort that list in descending order and get the indices of the corresponding images in the range [0, FSR * N)
	#	3. Add these images together. 
	if li_method == 'peak_pixel' and fsr < 1:
		sorted_idx = np.argsort(peak_pixel_vals)[::-1]	# Array holding indices of images
		N = np.ceil(fsr * N)
		image_stacked = (image_ref + np.sum(images_shifted[sorted_idx[:N]], 0)) / (N + 1)
	else:
		# Now, stacking the images. Need to change N if FSR < 1.
		image_stacked = (image_ref + np.sum(images_shifted, 0)) / (N + 1)	

	return image_stacked, rel_shift_idxs

####################################################################################################
def alignmentError(in_idxs, out_idxs,
	verbose=False):
	"""
		Compute the alignment errors arising in the Lucky Imaging shifting-and-stacking process given an input array of tip and tilt coordinates applied to the input images and the coordinates of the shifts applied in the shifting-and-stacking process.

		The total number of errors, the mean error and an array containing each alignment error is returned.i 
	"""
	N = in_idxs.shape[0]
	errs = np.zeros((N))
	n_errs = 0
	thresh = 0.1	# Threshold for misalignment
		
	for k in range(N):
		errs[k] = np.sqrt(np.power(in_idxs[k,0] - out_idxs[k,0],2) + np.power(in_idxs[k,1] - out_idxs[k,1],2))
		if errs[k] > thresh:
			n_errs += 1
	avg_err = np.mean(errs)
			
	if verbose:
		print('------------------------------------------------')
		print('Tip/tilt coordinates\nInput\t\tOutput\t\tError')
		print('------------------------------------------------')
		for k in range(N):
			print('(%6.2f,%6.2f)\t(%6.2f,%6.2f)\t%4.2f' % (in_idxs[k,0],in_idxs[k,1],out_idxs[k,0],out_idxs[k,1],errs[k]))
		print('------------------------------------------------')
		print('\t\t\tMean\t%4.2f' % avg_err)

	
	return n_errs, errs, avg_err

####################################################################################################
def _li_error_check(images, 
	image_ref = None,
	N = None):
	"""
		A private method to be used to check the inputs to the Lucky Imaging methods. 
	"""
	# Need to convert to float if necessary.
	if type(images.flatten()[0]) != np.float64:
		images = images.astype(np.float64)
	if image_ref is not None and type(image_ref.flatten()[0]) != np.float64:
		image_ref = image_ref.astype(np.float64)

	# Checking image dimensions.
	if len(images.shape) > 4:
		print("WARNING: for now, please only input a 3D array of images to shift and stack! I'm only going to operate on the first set of images...")
		images = np.squeeze(images[0])	
	
	if len(images.shape) == 3:
		if N and N > images.shape[0]:
			print("ERROR: if specified, N must be equal to or less than the length of the first dimension of the images array.")
			raise UserWarning
		if image_ref is None:
			# If the reference image is not specified, we use the first image in the array as the reference: 
			# i.e. we align all other images to images[0].
			if not N:
				N = images.shape[0]-1
			image_ref = np.copy(images[0])
			images = np.copy(images[1:])	# Only need to go through images 1:N-1.
		else:
			if image_ref.shape != images[0].shape:
				print("ERROR: if specified, reference image shape must be equal to input image stack shape.")
				raise UserWarning
			if not N:
				N = images.shape[0]			
	else:
		# Error: cannot shift and stack a single image!
		print("ERROR: cannot shift and stack a single image! Input array must have N > 1.")
		raise UserWarning

	return images, image_ref, N

####################################################################################################
def _centroid(image):
	""" Returns the centroid coordinates of an image. """
	height = image.shape[0]
	width = image.shape[1]
	x = np.arange(height)
	y = np.arange(width)
	X, Y = np.meshgrid(y,x)
	M_10 = np.sum((X * image).flatten())
	M_01 = np.sum((Y * image).flatten())
	M_00 = np.sum(image.flatten())

	centroid = np.asarray([M_01 / M_00, M_10 / M_00])

	return centroid