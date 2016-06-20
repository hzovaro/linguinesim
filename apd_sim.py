############################################################################################
#
# 	File:		apd_sim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		20/06/2016
#
#	Description:
#	Simulate images of galaxies imaged using the APD detector within the ANU 2.3 m telescope.
#
#	TO DO:
#	- find plate scales for HST images
#	- tip and tilt: convert the input sigma values to as
#	- implement poisson noise (not Gaussian)
#
############################################################################################

# import apdParameters as detector
# import anu23mParameters as telescope
# import cryoParameters as cryo
from apd_etc import *
import pdb

import numpy as np
import numpy.fft
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.signal
import scipy.misc
import matplotlib
matplotlib.rc('image', interpolation='none', cmap = 'gray')
plt.close('all')

#########################################################################################################
def getRawImages(fname):
	" Return an array of the image(s) stored in the FITS file fname. "
	hdulist = fits.open(fname)
	images_raw = hdulist[0].data
	hdulist.close()
	return images_raw	

#########################################################################################################
def resizeImagesToDetector(images_raw, source_plate_scale_as, dest_detector_size_px, dest_plate_scale_as):
	" Resize the images stored in array images_raw to the detector. "
	# 1. Get the original size and shape of the input images.
	if len(images_raw.shape) > 2:
		N = images_raw.shape[0]
		source_height_px = images_raw.shape[1]
		source_width_px = images_raw.shape[2]
	else:
		N = 1
		source_height_px = images_raw.shape[0]
		source_width_px = images_raw.shape[1]

	source_width_as = source_width_px * source_plate_scale_as
	source_height_as = source_height_px * source_plate_scale_as

	# 2. Resize the image to replicate the pixel scale of the detector:
	#	a. Get the angular extent of the source image:
	#		size(pixels on our detector) = size(of source, in as) / plate scale
	detector_height_px = dest_detector_size_px[0]
	detector_width_px = dest_detector_size_px[1]
	dest_width_px = source_width_as / dest_plate_scale_as
	dest_height_px = source_height_as / dest_plate_scale_as

	#	b. Rescale image to the appropriate size for our detector
	images = np.copy(images_raw)
	# images = scipy.misc.imresize(images, float(dest_width_px)/float(source_width_px))
	images = scipy.misc.imresize(images, (int(np.ceil(dest_height_px)), int(np.ceil(dest_width_px))))
	images = np.swapaxes(images,0,2)
	images = np.swapaxes(images,1,2)
	images = images.astype(np.float32)	# convert back to single precision floating point format


	#	c. Resize to detector.length * detector.width, maintaining the format. Padding if necessary
	if dest_height_px > detector_height_px:
		images = images[:, images.shape[1]/2-detector_height_px/2:images.shape[1]/2+detector_height_px/2+1, :]
		pad_height_top = 0
		pad_height_bottom = 0
	else:
		pad_height_top = np.floor((detector_height_px - images.shape[1])/2.).astype(np.int)
		pad_height_bottom = np.ceil((detector_height_px - images.shape[1])/2.).astype(np.int)

	if dest_width_px > detector.width_px:
		images = images[:, :, images.shape[2]/2-detector_width_px/2:images.shape[2]/2+detector_width_px/2+1]
		pad_width_left = 0
		pad_width_right = 0
	else: 
		pad_width_left = np.floor((detector_width_px - images.shape[2])/2.).astype(np.int)
		pad_width_right = np.ceil((detector_width_px - images.shape[2])/2.).astype(np.int)

	images = np.pad(images, ((0, 0), (pad_height_top, pad_height_bottom), (pad_width_left, pad_width_right)), mode='constant')

	return images	

# ##################################################################################
def getDiffractionLimitedImage(image_truth, wavelength, f_ratio, detector_size_px, l_px_m, 
	plotIt=False):
	" Convolve the PSF of a given telescope in a given band (J, H or K) with image_truth to simulate diffraction-limited imaging. "
	# Calculating the diffraction limit
	detector_height_px = detector_size_px[0]
	detector_width_px = detector_size_px[1]

	# Calculating the PSF
	y = np.arange(-detector_width_px/2, +detector_width_px/2, 1) * l_px_m
	x = np.arange(-detector_height_px/2, +detector_height_px/2, 1) * l_px_m
	X, Y = np.meshgrid(x, y)
	# x in units of m
	r = np.pi / wavelength / f_ratio * np.sqrt(np.power(X,2) + np.power(Y,2))
	# Calculating the PSF
	I_0 = 1
	psf = I_0 * np.power(2 * scipy.special.jv(1, r) / r, 2)
	psf[psf.shape[0]/2,psf.shape[1]/2] = I_0
	psf = np.swapaxes(psf,0,1)
	psf = psf.astype(np.float32)

	# Convolving the PSF and the truth image to obtain the simulated diffraction-limited image
	image_difflim = scipy.signal.fftconvolve(image_truth, psf, mode='same')

	if plotIt:
		plt.figure()

		plt.subplot(1,3,1)
		plt.imshow(psf)
		plt.title('Diffraction-limited PSF of telescope')
		plt.xlabel('x (pixels)')
		plt.ylabel('y (pixels)')

		plt.subplot(1,3,2)
		plt.imshow(image_truth)
		plt.title('Truth image')
		plt.xlabel('x (pixels)')
		plt.ylabel('y (pixels)')

		plt.subplot(1,3,3)
		plt.imshow(image_difflim)
		plt.title('Diffraction-limited image')
		plt.xlabel('x (pixels)')
		plt.ylabel('y (pixels)')

		plt.show()

	return image_difflim

#########################################################################################################
def addTurbulence(image, N, sigma_tt_px):
	" Add turbulence to an input `truth' image. Returns N copies of the input image with randomised turbulence added. "
	# Tip and tilt for now	
	# TODO add rotation
	# Shift by some random amount
	images = np.ndarray((N, image.shape[0], image.shape[1]))
	tt_idxs = np.ndarray((N, 2))
	for k in range(N):
		shift_width = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
		shift_height = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
		images[k] = np.roll(np.roll(image, shift_width, 0), shift_height, 1)
		tt_idxs[k] = [shift_width, shift_height]

	return (images, tt_idxs)

#########################################################################################################
"""
	ADDING NOISE:
		The Sigmga * t_exp that get output by the ETC are *what we should theoretically expect to measure 
		in every single pixel*. 
		However, because we are dealing with low numbers of photons here, we assume that the standard deviation goes as 
		the square root of the total numbers of pixels (Poisson statistics)
		Therefore the noise that we add should be a Gaussian *Centered* at Sigma * t_exp (the theoretical count rate)
		+/- sqrt(Sigma * t_exp) 
		The actual statistical properties of the nosie (i.e. the standard deviation) are actually determined by Poisson
		statistics!
"""
def addNoise(images,band,t_exp):
	" Add noise to an array of input images assuming an exposure time t_exp."
	# Creating an array in which to store the noisy images
	if len(images.shape) > 2:
		N = images.shape[0]
		height = images.shape[1]
		width = images.shape[2]
	else:
		N = 1
		height = images.shape[0]
		width = images.shape[1]

	noisyImages = np.copy(images)

	# Getting noise parameters from the ETC.
	etc_output = exposureTimeCalc(band,t_exp)

	# Adding noise to each image.
	for k in range(N):
		frame_sky = np.ones((height, width)) * etc_output['N_sky'] + np.random.randn(height, width) * etc_output['sigma_sky']
		frame_dark = np.ones((height, width)) * etc_output['N_dark'] + np.random.randn(height, width) * etc_output['sigma_dark']
		frame_cryo = np.ones((height, width)) * etc_output['N_cryo'] + np.random.randn(height, width) * etc_output['sigma_cryo']
		frame_RN = np.ones((height, width)) * etc_output['N_RN'] + np.random.randn(height, width) * etc_output['sigma_RN']
		if N > 1:
			noisyImages[k] += frame_sky + frame_cryo + frame_RN + frame_dark
		else:
			noisyImages += frame_sky + frame_cryo + frame_RN + frame_dark

	return noisyImages


#########################################################################################################
def shiftAndStack(images, 
	image_ref=None, 
	plotIt=False):
	" Shift and stack N images given in the array of N images."
	if (len(images.shape) > 2):
		if image_ref == None:
			# If the reference image is not specified, we use the first image in the array as the reference: 
			# i.e. we align all other images to images[0].
			N = images.shape[0]-1
			images = np.copy(images[1:])	# Only need to go through images 1:N-1.
			image_ref = np.copy(images[0])
		else:
			N = images.shape[0]			
	else:
		# Error: cannot shift and stack a single image!
		print 'ERROR: cannot shift and stack a single image! Input array must have N > 1.'
		return -1
	
	image_stacked = np.copy(image_ref)		# shifted-and-stacked image

	width = images[0].shape[0]
	height = images[0].shape[1]
	
	corrs = np.ndarray((N, 2 * width - 1, 2 * height - 1))	# Array to hold x-correlation results
	corr_peak_idxs = np.ndarray((N, 2))		# indices of the peak value in the x-correlation
	img_peak_idxs = np.ndarray((N, 2))		# shift in x and y computed from the x-correlation

	for k in range(N):
		# Cross-correlate image k with the reference image to find the tip and tilt.
		corrs[k] = scipy.signal.fftconvolve(image_ref, images[k][::-1,::-1])
		corr_peak_idxs[k] = np.unravel_index(np.argmax(corrs[k]), (2 * width - 1, 2 * height - 1))
		img_peak_idxs[k][0] = - corr_peak_idxs[k][0] + (width - 1)
		img_peak_idxs[k][1] = - corr_peak_idxs[k][1] + (height - 1)

		# Shift-and-stack the images.
		image_stacked += np.roll(np.roll(images[k], -img_peak_idxs[k][0].astype(int), 0), -img_peak_idxs[k][1].astype(int), 1)	

		# Plotting
		if plotIt:
			if k == 0:
				plt.figure()
				plt.subplot(1,3,1)
				plt.imshow(images[0],origin='lower')
				plt.subplot(1,3,2)
				scat2 = plt.scatter(0.0,0.0,c='r',s=20)
				plt.subplot(1,3,3)
				scat3 = plt.scatter(0.0,0.0,c='g',s=40)

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

	return image_stacked