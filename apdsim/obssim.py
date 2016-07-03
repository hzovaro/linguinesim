#########################################################################################################
#
# 	File:		obssim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A module for simulating imaging of objects using a given telescope and detector system.
#
#	Copyright (C) 2016 Anna Zovaro
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
#########################################################################################################
from __future__ import division
from apdsim import *

def resizeImagesToDetector(images_raw, source_plate_scale_as, dest_detector_size_px, dest_plate_scale_as,
	plotIt=False):
	" Resize the images stored in array images_raw with a given plate scale to a detector with given dimensions and plate scale. "
	print "Resizing image(s) to detector..."

	# 1. Get the original size and shape of the input images.
	images_raw, N, source_height_px, source_width_px = getImageSize(images_raw)
	source_width_as = source_width_px * source_plate_scale_as
	source_height_as = source_height_px * source_plate_scale_as

	# Getting the angular extent of the source image:
	# 	size(pixels on our detector) = size(of source, in as) / plate scale
	detector_height_px = dest_detector_size_px[0]
	detector_width_px = dest_detector_size_px[1]
	dest_width_px = source_width_as / dest_plate_scale_as
	dest_height_px = source_height_as / dest_plate_scale_as

	# Rescaling images to the appropriate size for our detector.
	images = np.ndarray((N, int(np.ceil(dest_height_px)), int(np.ceil(dest_width_px))))
	for k in range(N):
		im = Image.fromarray(images_raw[k])
		# NOTE: due to the way the Image package works, height and width indices are swapped
		im = im.resize((int(np.ceil(dest_width_px)), int(np.ceil(dest_height_px))), resample=PIL.Image.LANCZOS)
		images[k] = imageToArray(im)

	height_idx = 1	# Array index corresponding to image height.
	width_idx = 2	# Array index corresponding to image width.
		
	# Resizing to the size of the detector.
	if dest_height_px > detector_height_px:
		images = images[:, images.shape[height_idx]/2-detector_height_px/2:images.shape[height_idx]/2+detector_height_px/2, :]
		pad_height_top = 0
		pad_height_bottom = 0
	else:
		pad_height_top = np.floor((detector_height_px - images.shape[height_idx])/2.).astype(np.int)
		pad_height_bottom = np.ceil((detector_height_px - images.shape[height_idx])/2.).astype(np.int)

	if dest_width_px > detector_width_px:
		images = images[:, :, images.shape[width_idx]/2-detector_width_px/2:images.shape[width_idx]/2+detector_width_px/2]
		pad_width_left = 0
		pad_width_right = 0
	else: 
		pad_width_left = np.floor((detector_width_px - images.shape[width_idx])/2.).astype(np.int)
		pad_width_right = np.ceil((detector_width_px - images.shape[width_idx])/2.).astype(np.int)

	# Padding the resized images if necessary.
	images = np.pad(images, ((0, 0), (pad_height_top, pad_height_bottom), (pad_width_left, pad_width_right)), mode='constant')

	if plotIt:
		plt.figure()
		plt.subplot(1,2,1)
		plt.imshow(images_raw[0])
		plt.title('Input image')
		plt.subplot(1,2,2)
		plt.imshow(images[0])
		plt.title('Resized image')
		plt.suptitle('Resizing truth image to detector')
		plt.show()

	return np.squeeze(images)

# ##################################################################################
def getDiffractionLimitedImage(image_truth, wavelength, f_ratio, l_px_m, 
	detector_size_px=None,
	plotIt=False):
	" Convolve the PSF of a given telescope in a given band (J, H or K) with image_truth to simulate diffraction-limited imaging. "
	" It is assumed that the truth image has the appropriate plate scale of, but may be larger than, the detector. "
	" If the detector size is not given, then it is assumed that the input image and detector have the same dimensions. "
	print "Diffraction-limiting truth image(s)..."

	image_truth, N, height, width = getImageSize(image_truth)
	if detector_size_px != None:
		detector_height_px, detector_width_px = detector_size_px[0:2]
		if height < detector_height_px or width < detector_width_px:
			print "ERROR: the truth image must be larger than or the same size as the detector!"
			return -1		
	else:
		detector_height_px = height
		detector_width_px = width	

	# Calculating the PSF
	x = np.arange(-detector_height_px/2, +detector_height_px/2 + detector_height_px%2, 1) * l_px_m
	y = np.arange(-detector_width_px/2, +detector_width_px/2 + detector_width_px%2, 1) * l_px_m
	X, Y = np.meshgrid(x, y)
	# x in units of m
	r = np.pi / wavelength / f_ratio * np.sqrt(np.power(X,2) + np.power(Y,2))
	# Calculating the PSF
	I_0 = 1
	psf = I_0 * np.power(2 * special.jv(1, r) / r, 2)
	psf[psf.shape[0]/2,psf.shape[1]/2] = I_0	# removing the NaN in the centre of the image
	psf = np.swapaxes(psf,0,1)
	psf = psf.astype(np.float32)

	# Convolving the PSF and the truth image to obtain the simulated diffraction-limited image
	image_difflim = np.ndarray((N, height, width))
	for k in range(N):
		image_difflim[k] = signal.fftconvolve(image_truth[k], psf, mode='same')

	if plotIt:
		plt.figure()
		plt.subplot(1,3,1)
		plt.imshow(psf)
		plt.title('Diffraction-limited PSF of telescope')
		plt.subplot(1,3,2)
		plt.imshow(image_truth[0])
		plt.title('Truth image')
		plt.subplot(1,3,3)
		plt.imshow(image_difflim[0])
		plt.title('Diffraction-limited image')
		plt.suptitle('Diffraction-limiting image')
		plt.show()

	return np.squeeze(image_difflim)

#########################################################################################################
def getSeeingLimitedImage(images, seeing_diameter_as, plate_scale_as,
	padFactor=1,
	plotIt=False):
	" Convolve a Gaussian PSF with an input image to simulate seeing with a FWHM of seeing_diameter_as. "
	print "Seeing-limiting image(s)..."

	images, N, height, width = getImageSize(images)

	# Padding the source image.
	pad_ud = height / padFactor / 2
	pad_lr = width / padFactor / 2
	
	# If the image dimensions are odd, need to ad an extra row/column of zeros.
	# image_padded = np.pad(images[0], ((pad_ud,pad_ud + height % 2),(pad_lr,pad_lr + width % 2)), mode='constant')
	# conv_height = image_padded.shape[0]
	# conv_width = image_padded.shape[1]
	conv_height = 2 * pad_ud + height + (height % 2)
	conv_width = 2 * pad_lr + width + (width % 2)

	# Generate a Gaussian kernel.
	kernel = np.zeros(image_padded.shape)
	y_as = np.arange(-conv_width /2, +conv_width/2, 1) * plate_scale_as
	x_as = np.arange(-conv_height/2, +conv_height/2, 1) * plate_scale_as
	X, Y = np.meshgrid(x_as, y_as)
	sigma = seeing_diameter_as / (2 * np.sqrt(2 * np.log(2)))
	kernel = np.exp(-(np.power(X, 2) + np.power(Y, 2)) / (2 * np.power(sigma,2)))
	kernel = np.pad(kernel, ((pad_ud, pad_ud + height % 2), (pad_lr, pad_lr + width % 2)), mode='constant')

	# Convolving the kernal with the image.
	image_seeing_limited = np.ndarray((N, conv_height, conv_width))
	image_seeing_limited_cropped = np.ndarray((N, height, width))

	for k in range(N):
		image_padded = np.pad(images[k], ((pad_ud,pad_ud + height % 2),(pad_lr,pad_lr + width % 2)), mode='constant')
		image_seeing_limited[k] = signal.fftconvolve(image_padded, kernel, mode='same')
		image_seeing_limited_cropped[k] = image_seeing_limited[pad_ud : height + pad_ud, pad_lr : width + pad_lr]		

	if plotIt:
		plt.figure()
		plt.subplot(2,2,1)
		plt.imshow(image)
		plt.title('Input image')
		plt.subplot(2,2,2)
		plt.imshow(kernel)
		plt.title('Kernel')
		plt.subplot(2,2,3)
		plt.imshow(image_seeing_limited[0])
		plt.title('Convolved image')
		plt.subplot(2,2,4)
		plt.imshow(image_seeing_limited_cropped[0])
		plt.title('Cropped, convolved image')

	return np.squeeze(image_seeing_limited_cropped)

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
	print ("Adding noise to image(s)...")

	# Creating an array in which to store the noisy images
	images, N, height, width = getImageSize(images)

	noisyImages = np.copy(images)

	# Getting noise parameters from the ETC.
	etc_output = exposureTimeCalc(band,t_exp)

	# Adding noise to each image.
	for k in range(N):
		frame_sky = np.ones((height, width)) * etc_output['N_sky'] + np.random.randn(height, width) * etc_output['sigma_sky']
		frame_dark = np.ones((height, width)) * etc_output['N_dark'] + np.random.randn(height, width) * etc_output['sigma_dark']
		frame_cryo = np.ones((height, width)) * etc_output['N_cryo'] + np.random.randn(height, width) * etc_output['sigma_cryo']
		frame_RN = np.ones((height, width)) * etc_output['N_RN'] + np.random.randn(height, width) * etc_output['sigma_RN']
		noisyImages[k] += frame_sky + frame_cryo + frame_RN + frame_dark		

	return (np.squeeze(noisyImages), etc_output)
