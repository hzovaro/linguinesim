############################################################################################
#
# 	File:		apd_sim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		27/06/2016
#
#	Description:
#	Simulate images of galaxies imaged using the APD detector within the ANU 2.3 m telescope.
#
#	TO DO:
#	- fix bug in resizing image?
#	- find plate scales for HST images
#	- tip and tilt: convert the input sigma values to as
#	- implement poisson noise (not Gaussian)
#
############################################################################################

from apd_etc import *
import pdb
from scipy.ndimage.interpolation import shift
import PIL
from PIL import Image
import numpy as np
import numpy.fft
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.signal
import scipy.misc
import matplotlib
from matplotlib.colors import LogNorm 

matplotlib.rc('image', interpolation='none', cmap = 'binary')
plt.close('all')

#########################################################################################################
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

	if plotIt:
		if len(images_raw.shape) > 2:
			N = images_raw.shape[0]
			for k in range(N):
				plt.figure()
				plt.imshow(images_raw[k])
				plt.title('Raw image %d from FITS file' % k + 1)
		else:
			plt.figure()
			plt.imshow(images_raw)
			plt.title('Raw image from FITS file')
		plt.show()

	return images_raw, hdulist	

#########################################################################################################
def resizeImagesToDetector(images_raw, source_plate_scale_as, dest_detector_size_px, dest_plate_scale_as,
	plotIt=False):
	" Resize the images stored in array images_raw with a given plate scale to a detector with given dimensions and plate scale. "
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

	# Getting the angular extent of the source image:
	# 	size(pixels on our detector) = size(of source, in as) / plate scale
	detector_height_px = dest_detector_size_px[0]
	detector_width_px = dest_detector_size_px[1]
	dest_width_px = source_width_as / dest_plate_scale_as
	dest_height_px = source_height_as / dest_plate_scale_as

	# Rescaling images to the appropriate size for our detector.
	if N > 1:
		images = np.ndarray([N, int(np.ceil(dest_height_px)), int(np.ceil(dest_width_px))])
		for k in range(N):
			im = Image.fromarray(images_raw[k])
			# NOTE: due to the way the Image package works, height and width indices are swapped
			im = im.resize((int(np.ceil(dest_width_px)), int(np.ceil(dest_height_px))), resample=PIL.Image.LANCZOS)
			images[k] = imageToArray(im)
	else:
		images = np.ndarray([int(np.ceil(dest_height_px)), int(np.ceil(dest_width_px))])
		im = Image.fromarray(images_raw)
		# NOTE: due to the way the Image package works, height and width indices are swapped
		im = im.resize((int(np.ceil(dest_width_px)), int(np.ceil(dest_height_px))), resample=PIL.Image.LANCZOS)
		images = imageToArray(im)

	if N > 1:
		height_idx = 1	# Array index corresponding to image height.
		width_idx = 2	# Array index corresponding to image width.
		# Swapping axes around because imresize does weird things.
		images = np.swapaxes(images,0,2)
		images = np.swapaxes(images,1,2)
	else:
		height_idx = 0	# Array index corresponding to image height.
		width_idx = 1	# Array index corresponding to image width.
		# Swapping axes around because imresize does weird things.
		# images = np.fliplr(images)
		# images = np.swapaxes(images,0,1)
	
	# Resizing to the size of the detector.
	if dest_height_px > detector_height_px:
		if N > 1:
			images = images[:, images.shape[height_idx]/2-detector_height_px/2:images.shape[height_idx]/2+detector_height_px/2, :]
		else:
			images = images[images.shape[height_idx]/2-detector_height_px/2:images.shape[height_idx]/2+detector_height_px/2, :]
		pad_height_top = 0
		pad_height_bottom = 0
	else:
		pad_height_top = np.floor((detector_height_px - images.shape[height_idx])/2.).astype(np.int)
		pad_height_bottom = np.ceil((detector_height_px - images.shape[height_idx])/2.).astype(np.int)

	if dest_width_px > detector_width_px:
		if N > 1:
			images = images[:, :, images.shape[width_idx]/2-detector_width_px/2:images.shape[width_idx]/2+detector_width_px/2]
		else:
			images = images[:, images.shape[width_idx]/2-detector_width_px/2:images.shape[width_idx]/2+detector_width_px/2]
		pad_width_left = 0
		pad_width_right = 0
	else: 
		pad_width_left = np.floor((detector_width_px - images.shape[width_idx])/2.).astype(np.int)
		pad_width_right = np.ceil((detector_width_px - images.shape[width_idx])/2.).astype(np.int)

	# Padding the resized images if necessary.
	if N > 1:
		images = np.pad(images, ((0, 0), (pad_height_top, pad_height_bottom), (pad_width_left, pad_width_right)), mode='constant')
	else:
		images = np.pad(images, ((pad_height_top, pad_height_bottom), (pad_width_left, pad_width_right)), mode='constant')

	if plotIt:
		plt.figure()
		if N == 1:
			plt.subplot(1,2,1)
			plt.imshow(images_raw)
			plt.title('Input image')
			plt.subplot(1,2,2)
			plt.imshow(images)
			plt.title('Resized image')
		else:
			plt.subplot(1,2,1)
			plt.imshow(images_raw[0])
			plt.title('Input image')
			plt.subplot(1,2,2)
			plt.imshow(images[0])
			plt.title('Resized image')
		plt.suptitle('Resizing truth image to detector')
		plt.show()

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
		plt.subplot(1,3,2)
		plt.imshow(image_truth)
		plt.title('Truth image')
		plt.subplot(1,3,3)
		plt.imshow(image_difflim)
		plt.title('Diffraction-limited image')
		plt.suptitle('Diffraction-limiting image')
		plt.show()

	return image_difflim

#########################################################################################################
def getSeeingLimitedImage(image, seeing_diameter_as, plate_scale_as,
	padFactor=1,
	plotIt=False):
	detector_height_px = image.shape[0]
	detector_width_px = image.shape[1]

	# Padding the source image.
	pad_ud = detector_height_px / padFactor / 2
	pad_lr = detector_width_px / padFactor / 2
	
	# If the image dimensions are odd, need to ad an extra row/column of zeros.
	image_padded = np.pad(image, ((pad_ud,pad_ud + detector_height_px % 2),(pad_lr,pad_lr + detector_width_px % 2)), mode='constant')
	conv_height = image_padded.shape[0]
	conv_width = image_padded.shape[1]

	# Generate a Gaussian kernel.
	kernel = np.zeros(image_padded.shape)
	y_as = np.arange(-conv_width /2, +conv_width/2, 1) * plate_scale_as
	x_as = np.arange(-conv_height/2, +conv_height/2, 1) * plate_scale_as
	X, Y = np.meshgrid(x_as, y_as)
	sigma = seeing_diameter_as / (2 * np.sqrt(2 * np.log(2)))
	kernel = np.exp(-(np.power(X, 2) + np.power(Y, 2)) / (2 * np.power(sigma,2)))
	kernel = np.pad(kernel, ((pad_ud, pad_ud + detector_height_px % 2), (pad_lr, pad_lr + detector_width_px % 2)), mode='constant')

	# Convolving the kernal with the image.
	image_seeing_limited = scipy.signal.fftconvolve(image_padded, kernel, mode='same')
	image_seeing_limited_cropped = image_seeing_limited[pad_ud : detector_height_px + pad_ud, pad_lr : detector_width_px + pad_lr]

	if plotIt:
		plt.figure()
		plt.subplot(2,2,1)
		plt.imshow(image)
		plt.title('Input image')
		plt.subplot(2,2,2)
		plt.imshow(kernel)
		plt.title('Kernel')
		plt.subplot(2,2,3)
		plt.imshow(image_seeing_limited)
		plt.title('Convolved image')
		plt.subplot(2,2,4)
		plt.imshow(np.log(image_seeing_limited_cropped))
		plt.title('Cropped, convolved image')

	return image_seeing_limited_cropped

#########################################################################################################
def addTurbulence(image, N, sigma_tt_px,
	crop_tt=None):
	" Add turbulence to an input `truth' image. Returns N copies of the input image with randomised turbulence added. "

	# Tip and tilt for now	
	height = image.shape[0]
	width = image.shape[1]

	# Output array of images
	if crop_tt == None:
		images_tt = np.ndarray((N, height, width))
	else:
		if type(crop_tt) == int:
			images_tt = np.ndarray((N, height - 2 * crop_tt, width - 2 * crop_tt))	
		else:
			images_tt = np.ndarray((N, height - 2 * crop_tt[0], width - 2 * crop_tt[1]))

	# Array to hold the tip/tilt offsets
	tt_idxs = np.ndarray((N, 2))
	
	for k in range(N):
		shift_height = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
		shift_width = np.ceil(np.random.randn() * sigma_tt_px).astype(int)
		image_tt = shift(image, (shift_height, shift_width))
		tt_idxs[k] = [shift_height, shift_width]

		# Cropping the image if necessary
		# (left, upper, right, lower)
		if crop_tt == None:
			images_tt[k] = image_tt
		else:
			# if type(crop_tt) == int:
				# cropIdx = (crop_tt, crop_tt, - crop_tt + width, - crop_tt + height)
			# else:
				# cropIdx = (crop_tt[1], crop_tt[0], - crop_tt[1] + width, - crop_tt[0] + height)
			# images_tt[k] = rotateAndCrop(image_tt, angle=0., cropArg=cropIdx)
			images_tt[k] = rotateAndCrop(image_tt, angle=0., cropArg=crop_tt)

	return (images_tt, tt_idxs)

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

	return (noisyImages, etc_output)


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

	height = images[0].shape[0]
	width = images[0].shape[1]
	
	corrs = np.ndarray((N, 2 * height - 1, 2 * width - 1))	# Array to hold x-correlation results
	corr_peak_idxs = np.ndarray((N, 2))		# indices of the peak value in the x-correlation
	img_peak_idxs = np.ndarray((N, 2))		# shift in x and y computed from the x-correlation

	for k in range(N):
		# Cross-correlate image k with the reference image to find the tip and tilt.
		corrs[k] = scipy.signal.fftconvolve(image_ref, images[k][::-1,::-1])
		corr_peak_idxs[k] = np.unravel_index(np.argmax(corrs[k]), (2 * height - 1, 2 * width - 1))
		img_peak_idxs[k][0] = - corr_peak_idxs[k][0] + (height - 1)
		img_peak_idxs[k][1] = - corr_peak_idxs[k][1] + (width - 1)

		# Shift-and-stack the images.
		# image_stacked += np.roll(np.roll(images[k], -img_peak_idxs[k][0].astype(int), 0), -img_peak_idxs[k][1].astype(int), 1)	
		image_stacked += shift(images[k], (-img_peak_idxs[k][0].astype(int), -img_peak_idxs[k][1].astype(int)))

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
	" Rotate and crop an array of N images stored in ndarray image_in_array counterclockwise by a given 	angle and then crop the image using coordinates (left, upper, right, lower) "
	if len(image_in_array.shape) > 2:
		N = image_in_array.shape[0]	
		height = image_in_array.shape[1]
		width = image_in_array.shape[2]	
	else:
		N = 1
		height = image_in_array.shape[0]
		width = image_in_array.shape[1]

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
	if N == 1:
		image = Image.fromarray(image_in_array)
		# Rotate.
		if angle != 0.:
			image = image.rotate(angle)
		# Crop.
		image = image.crop(cropIdx)
		# Convert back to an array.
		image_out_array = imageToArray(image)
	else:
		image_out_array = np.ndarray((N, height, width))
		for k in range(N):
			image = Image.fromarray(image_in_array[k])
			# Rotate.
			if angle != 0.:
				image = image.rotate(angle)
			# Crop.
			image = image.crop(cropIdx)
			# Convert back to an array.
			image_out_array[k] = imageToArray(image)

	if plotIt:
		plt.figure()
		if N > 1:			
			plt.subplot(1,2,1)
			plt.imshow(image_in_array[0])
			plt.title('Input image')
			plt.subplot(1,2,2)
			plt.imshow(image_out_array[0])
			plt.title('Output image')
		else:
			plt.subplot(1,2,1)
			plt.imshow(image_in_array)
			plt.title('Input image')
			plt.subplot(1,2,2)
			plt.imshow(image_out_array)
			plt.title('Output image')
		plt.suptitle('Rotating and cropping image')
		plt.show()

	return image_out_array

