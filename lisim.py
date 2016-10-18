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
from matplotlib.cbook import is_numlike
from matplotlib.mlab import normpdf
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
def luckyImage(im, psf, scale_factor, t_exp, final_sz,
	tt = np.array([0, 0]),
	im_star = None,					# Image of a star at the pre-downsizing plate scale to add to the image.					
	noise_frame_pre_gain = 0,		# Noise added before gain multiplication, e.g. sky background, emission from telescope, etc. Must have shape final_sz.
	noise_frame_post_gain = 0,		# Noise added after gain multiplication, e.g. read noise. Must have shape final_sz.
	gain = 1,						# Detector gain.
	plate_scale_as_px_conv = 1,		# Only used for plotting.
	plate_scale_as_px = 1,			# Only used for plotting.
	plotIt=False):
	""" 
		This function can be used to generate a short-exposure 'lucky' image that can be input to the Lucky Imaging algorithms.
			Input: 	one 'raw' countrate image of a galaxy; one PSF with which to convolve it (at the same plate scale)
			Output: a 'Lucky' exposure. 			
			Process: convolve with PSF --> resize to detector --> add tip and tilt (from a premade vector of tip/tilt values) --> convert to counts --> add noise --> subtract the master sky/dark current. 
	"""	
	# Convolve with PSF.
	im_raw = im
	im_convolved = obssim.convolvePSF(im_raw, psf)
		
	# Add a star to the field. We need to add the star at the convolution plate scale BEFORE we resize down because of the tip-tilt adding step!
	if is_numlike(im_star):
		if im_star.shape != im_convolved.shape:
			print("ERROR: the input image of the star MUST have the same size and plate scale as the image of the galaxy after convolution!")
			raise UserWarning
		im_convolved += im_star

	# Resize to detector (+ edge buffer).
	im_resized = obssim.resizeImageToDetector(image_raw = im_convolved, source_plate_scale_as = 1, dest_plate_scale_as = scale_factor, conserve_pixel_sum=True)	
	# Add tip and tilt. To avoid edge effects, max(tt) should be less than or equal to the edge buffer.
	edge_buffer_px = (im.shape[0] - final_sz[0]) / 2
	if edge_buffer_px > 0 and max(tt) > edge_buffer_px:
		print("WARNING: the edge buffer is less than the supplied tip and tilt by a margin of {:.2f} pixels! Shifted image will be clipped.".format(np.abs(edge_buffer_px - max(tt))))
	im_tt = obssim.addTipTilt_single(image = im_resized, tt_idxs = tt)[0]	
	# Crop back down to the detector size.
	if edge_buffer_px > 0:
		im_tt = imutils.centreCrop(im_tt, final_sz)	
	# Convert to counts.
	im_counts = etcutils.expectedCount2count(im_tt, t_exp = t_exp)
	# Add the pre-gain noise. 
	im_noisy_pre_gain = im_counts + noise_frame_pre_gain
	# Multiply by the detector gain.
	im_noisy = im_noisy_pre_gain * gain 
	# Add the post-gain noise.
	im_noisy += noise_frame_post_gain

	if plotIt:
		# Plotting
		mu.newfigure(1,4)
		plt.suptitle('Convolving input image with PSF and resizing to detector')
		mu.astroimshow(im=im_raw, title='Raw input image (electrons/s)', plate_scale_as_px = plate_scale_as_px_conv, colorbar_on=True, subplot=141)
		mu.astroimshow(im=psf, title='Point spread function (normalised)', plate_scale_as_px = plate_scale_as_px_conv, colorbar_on=True, subplot=142)
		mu.astroimshow(im=im_convolved, title='Convolved with PSF (electrons/s)', plate_scale_as_px = plate_scale_as_px_conv, colorbar_on=True, subplot=143)
		mu.astroimshow(im=im_resized, title='Resized to detector plate scale (electrons/s)', plate_scale_as_px=plate_scale_as_px, colorbar_on=True, subplot=144)

		mu.newfigure(1,4)
		plt.suptitle('Adding tip and tilt, converting to integer counts and adding noise')		
		mu.astroimshow(im=im_tt, title='Atmospheric tip and tilt added (electrons/s)', plate_scale_as_px=plate_scale_as_px, colorbar_on=True, subplot=141)
		mu.astroimshow(im=im_counts, title=r'Cropped, converted to integer counts (electrons)', plate_scale_as_px=plate_scale_as_px, colorbar_on=True, subplot=142)
		mu.astroimshow(im=im_noisy, title='Noise added (electrons)', plate_scale_as_px=plate_scale_as_px, colorbar_on=True, subplot=143)
		plt.subplot(1,4,4)
		x = np.linspace(-im_tt.shape[0]/2, +im_tt.shape[0]/2, im_tt.shape[0]) * plate_scale_as_px
		plt.plot(x, im_tt[:, im_tt.shape[1]/2], 'g', label='Electron count rate')
		plt.plot(x, im_counts[:, im_tt.shape[1]/2], 'b', label='Converted to integer counts ($t_{exp} = %.2f$ s)' % t_exp)
		plt.plot(x, im_noisy[:, im_tt.shape[1]/2], 'r', label='Noise added')
		plt.xlabel('arcsec')
		plt.ylabel('Pixel value (electrons)')
		plt.title('Linear profiles')
		plt.axis('tight')
		plt.legend()
		plt.show()

	return im_noisy

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

	# # Thresholding the image
	# image_subtracted_bg = np.copy(image)
	# image_subtracted_bg[image<1.5*np.mean(image.flatten())] = 0

	img_peak_idx = _centroid(image)

	# Shift the image by the relative amount.
	rel_shift_idx = (img_ref_peak_idx - img_peak_idx)
	image_shifted = scipy.ndimage.interpolation.shift(image, rel_shift_idx)

	return image_shifted, -rel_shift_idx

####################################################################################################
def shift_xcorr(image, image_ref, buff, subPixelShift):
	if type(image) == list:
		image = np.array(image)
	
	# Subtracting the mean of each image
	image_subtracted_bg = image - np.mean(image.flatten())
	image_ref_subtracted_bg = image_ref - np.mean(image_ref.flatten())

	height, width = image.shape
	if fftwconvolve.NTHREADS==0:
		corr = scipy.signal.fftconvolve(image_ref_subtracted_bg, image_subtracted_bg[::-1,::-1], 'same')
	else:
		corr = fftwconvolve.fftconvolve(image_ref_subtracted_bg, image_subtracted_bg[::-1,::-1], 'same')
	corr /= max(corr.flatten())	# The fitting here does not work if the pixels have large values!
	
	if subPixelShift: 
		# Fitting a Gaussian.
		Y, X = np.mgrid[-(height-2*buff)/2:(height-2*buff)/2, -(width-2*buff)/2:(width-2*buff)/2]
		x_peak, y_peak = np.unravel_index(np.argmax(corr), corr.shape)
		try:		
			p_init = astropy.modeling.models.Gaussian2D(x_mean=X[x_peak,y_peak],y_mean=Y[x_peak,y_peak],x_stddev=5.,y_stddev=5.,amplitude=np.max(corr.flatten()))
		except:			
			p_init = astropy.modeling.models.Gaussian2D(x_mean=x_peak,y_mean=y_peak,x_stddev=1.,y_stddev=1.,amplitude=1.)
		fit_p = astropy.modeling.fitting.LevMarLSQFitter()
		p_fit = fit_p(p_init, X, Y, corr[buff:height-buff, buff:width-buff])
		rel_shift_idx = (p_fit.y_mean.value, p_fit.x_mean.value)	# NOTE: the indices have to be swapped around here for some reason!		
	else:
		rel_shift_idx = np.unravel_index(np.argmax(corr), corr.shape)
		rel_shift_idx = (rel_shift_idx[0] - height/2, rel_shift_idx[1] - width/2)
	
	# mu.newfigure(1,5)
	# mu.astroimshow(im=image, title='Input image', subplot=151)
	# mu.astroimshow(im=image_ref, title='Reference image', subplot=152)
	# mu.astroimshow(im=corr, title='Cross-correlation', subplot=153)
	# mu.astroimshow(im=p_init(X,Y), title='p_init', subplot=154)
	# mu.astroimshow(im=p_fit(X,Y), title='p_fit', subplot=155)
	# plt.show()
	# pdb.set_trace()

	image_shifted = scipy.ndimage.interpolation.shift(image, rel_shift_idx)	

	return image_shifted, tuple(-x for x in rel_shift_idx)

####################################################################################################
def shift_gaussfit(image, img_ref_peak_idx):
	if type(image) == list:
		image = np.array(image)

	# Subtracting the mean of the input image
	image_subtracted_bg = image - np.mean(image.flatten())

	# Fitting a Gaussian to the mean-subtracted image.
	peak_idx = _gaussfit_peak(image_subtracted_bg)	
	rel_shift_idx = -(peak_idx - img_ref_peak_idx)

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
	stacking_method = 'average',
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
	li_method = li_method.lower()
	if li_method == 'cross-correlation':
		shift_fun = partial(shift_xcorr, image_ref=image_ref, buff=buff, subPixelShift=subPixelShift)	
	elif li_method == 'gaussian fit':
		img_ref_peak_idx = _gaussfit_peak(image_ref - np.mean(image_ref.flatten()))
		shift_fun = partial(shift_gaussfit, img_ref_peak_idx=img_ref_peak_idx)
	elif li_method == 'peak pixel':
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
		print("ERROR: invalid Lucky Imaging method '{}' specified; must be 'cross-correlation', 'peak pixel', 'centroid' or 'Gaussian fit' for now...".format(li_method))
		raise UserWarning

	# In here, want to parallelise the processing for *each image*. So make shift functions that work on a single image and return the shifted image, then stack it out here.
	
	tic = time.time()
	if mode == 'parallel':
		# Setting up to execute in parallel.
		images = images.tolist()	# Need to convert the image array to a list.

		# Executing in parallel.
		if fftwconvolve.NTHREADS == 0:
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
		if li_method == 'peak pixel' and fsr < 1:
			peak_pixel_vals = np.array(zip(*results)[2])

	elif mode == 'serial':
		# Loop through each image individually.
		images_shifted = np.zeros( (N, image_ref.shape[0], image_ref.shape[1]) )	
		rel_shift_idxs = np.zeros( (N, 2) )
		for k in range(N):
			if li_method == 'peak pixel':
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
	if li_method == 'peak pixel' and fsr < 1:
		sorted_idx = np.argsort(peak_pixel_vals)[::-1]	# Array holding indices of images
		N = np.ceil(fsr * N)
		# Is averaging the best way to do this? Probably not...
		if stacking_method == 'median combine':
			arr = np.ndarray((1, image_ref.shape[0], image_ref.shape[1]))
			arr[0,:] = image_ref
			image_ref = arr
			image_stacked = obssim.medianCombine(np.concatenate((image_ref, images_shifted[sorted_idx[:N]])))
		elif stacking_method == 'average':
			image_stacked = (image_ref + np.sum(images_shifted[sorted_idx[:N]], 0)) / (N + 1)
	else:
		# Now, stacking the images. Need to change N if FSR < 1.
		if stacking_method == 'median combine':			
			arr = np.ndarray((1, image_ref.shape[0], image_ref.shape[1]))
			arr[0,:] = image_ref
			image_ref = arr
			image_stacked = obssim.medianCombine(np.concatenate((image_ref, images_shifted)))
		elif stacking_method == 'average':
			image_stacked = (image_ref + np.sum(images_shifted, 0)) / (N + 1)	

	return image_stacked, rel_shift_idxs

####################################################################################################
def alignmentError(in_idxs, out_idxs, opticalsystem,
	li_method='',
	plotHist=True,
	verbose=True):
	"""
		Compute the alignment errors arising in the Lucky Imaging shifting-and-stacking process given an input array of tip and tilt coordinates applied to the input images and the coordinates of the shifts applied in the shifting-and-stacking process.
	"""
	N = in_idxs.shape[0]
	errs_as = np.zeros( (N, 3) )
		
	for k in range(N):
		errs_as[k, 0] = (in_idxs[k,0] - out_idxs[k,0]) * opticalsystem.plate_scale_as_px
		errs_as[k, 1] = (in_idxs[k,1] - out_idxs[k,1]) * opticalsystem.plate_scale_as_px
		errs_as[k, 2] = np.sqrt(errs_as[k, 0]**2 + errs_as[k, 1]**2)
	errs_px = errs_as / opticalsystem.plate_scale_as_px
			
	# Print the alignment errors to screen.
	if verbose:
		print('------------------------------------------------')
		print('Tip/tilt coordinates\nInput\t\tOutput\t\tError\tError (arcsec)')
		print('------------------------------------------------')
		for k in range(N):
			print('(%6.2f,%6.2f)\t(%6.2f,%6.2f)\t%4.2f\t%4.2f' % (in_idxs[k,0],in_idxs[k,1],out_idxs[k,0],out_idxs[k,1],errs_px[k,2],errs_as[k,2]))
		print('------------------------------------------------')
		print('\t\t\tMean\t%4.2f' % np.mean(errs_as))

	if plotHist:
		plotErrorHistogram(errs_as)
	
	return errs_as

####################################################################################################
def plotErrorHistogram(errs_as)
	x_errs_as = errs_as[:,0]
	y_errs_as = errs_as[:,1]
	# Plot a pretty shistogram showing the distribution of the alignment errors, and fit a Gaussian to them.
	range_as = 2 * max(max(np.abs(y_errs_as)), max(np.abs(x_errs_as)))
	nbins = int(N / 10)
	mu.newfigure(1.5,1)
	plt.suptitle('{} Lucky Imaging shifting-and-stacking alignment errors'.format(li_method))

	plt.subplot(211)
	plt.hist(x_errs_as, bins=nbins, range=(-range_as/2,+range_as/2), normed=True)
	mean_x = np.mean(x_errs_as)
	sigma_x = np.sqrt(np.var(x_errs_as))
	x = np.linspace(-range_as/2, range_as/2, 100)
	plt.plot(x, normpdf(x,mean_x,sigma_x), 'r', label=r'$\sigma_x$ = %.4f"' % (sigma_x))
	plt.title(r'$x$ alignment error')
	plt.xlabel('arcsec')
	plt.legend()

	plt.subplot(212)
	plt.hist(y_errs_as, bins=nbins, range=(-range_as/2,+range_as/2), normed=True)
	mean_y = np.mean(y_errs_as)
	sigma_y = np.sqrt(np.var(y_errs_as))
	y = np.linspace(-range_as/2, range_as/2, 100)
	plt.plot(y, normpdf(y,mean_y,sigma_y), 'r', label=r'$\sigma_y$ = %.4f"' % (sigma_y))
	plt.title(r'$y$ alignment error')
	plt.xlabel('arcsec')
	plt.legend()	

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

####################################################################################################
def _gaussfit_peak(image):
	""" Returns the coordinates of the peak of a 2D Gaussian fitted to the image. """
	height, width = image.shape

	# Fitting a Gaussian.
	Y, X = np.mgrid[-(height)/2:(height)/2, -(width)/2:(width)/2]
	try:		
		p_init = astropy.modeling.models.Gaussian2D(x_stddev=1.,y_stddev=1.)
	except:
		p_init = astropy.modeling.models.Gaussian2D(x_mean=1.,y_mean=1.,x_stddev=1.,y_stddev=1.,amplitue=1.)
	fit_p = astropy.modeling.fitting.LevMarLSQFitter()
	p_fit = fit_p(p_init, X, Y, image)		
	
	peak_idx = np.array([p_fit.y_mean.value, p_fit.x_mean.value])	# NOTE: the indices have to be swapped around here for some reason!		

	return peak_idx
