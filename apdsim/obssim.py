####################################################################################################
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
####################################################################################################
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

def airyDisc(wavelength_m, f_ratio, l_px_m, detector_size_px,
	x_offset=0,
	y_offset=0,
	P_0=1,
	plotIt=False):
	"""
		Returns the PSF of an optical system with a circular aperture given the f ratio, pixel and detector size at a given wavelength_m.

		If desired, an offset (measured from the top left corner of the detector) can be specified.

		The PSF is normalised such that the sum of every pixel in the PSF (extended to infinity) is equal to P_0 (unity by default), where P_0 is the total energy incident upon the telescope aperture. 

		P_0 represents the *ideal* total energy in the airy disc (that is, the total energy incident upon the telescope aperture), whilst P_sum measures the actual total energy in the image (i.e. the pixel values). 
	"""
	# Calculating the PSF
	detector_height_px, detector_width_px = detector_size_px[0:2]
	dx = detector_height_px/2 - x_offset
	dy = detector_width_px/2 - y_offset
	x = np.arange(-detector_height_px//2, +detector_height_px//2 + detector_height_px%2, 1) + dx
	y = np.arange(-detector_width_px//2, +detector_width_px//2 + detector_width_px%2, 1) + dy
	x *= l_px_m
	y *= l_px_m
	X, Y = np.meshgrid(x, y)
	# x in units of m
	r_prime = np.sqrt(np.power(X,2) + np.power(Y,2))
	r = np.pi / wavelength_m / f_ratio * r_prime
	# Central intensity (W m^-2)
	I_0 = P_0 * np.pi / 4 / wavelength_m / wavelength_m / f_ratio / f_ratio
	# Calculating the PSF
	airyDisc = np.power((2 * special.jv(1, r) / r), 2) * I_0 
	nan_idx = np.where(np.isnan(airyDisc))
	if nan_idx[0].shape != (0,):
		airyDisc[nan_idx[0][0],nan_idx[1][0]] = I_0 # removing the NaN in the centre of the image if necessary
	airyDisc *= l_px_m * l_px_m
	airyDisc = np.swapaxes(airyDisc,0,1)
	airyDisc = airyDisc.astype(np.float64)
	P_sum = sum(airyDisc.flatten())

	if plotIt:
		plt.figure(figsize=(FIGSIZE,FIGSIZE))
		# plt.imshow(psf, norm=LogNorm())
		plt.imshow(airyDisc)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Point spread function')
		plt.show()

	return airyDisc, P_0, P_sum, I_0

####################################################################################################
def psfKernel(wavelength_m, f_ratio, l_px_m, detector_size_px,
	plotIt=False):
	"""
		Returns an Airy disc PSF corresponding to an optical system with a given f ratio, pixel size and detector size at a specified wavelength_m.

		Note: if you are using this routine to generate a kernel with which to convolve a truth image, make sure the kernel is twice the size of the truth image!
	"""	
	return airyDisc(wavelength_m=wavelength_m, f_ratio=f_ratio, l_px_m=l_px_m, detector_size_px=detector_size_px, plotIt=plotIt)	

####################################################################################################
def resizeImagesToDetector(images_raw, source_plate_scale_as, dest_detector_size_px, dest_plate_scale_as,
	plotIt=False):
	" Resize the images stored in array images_raw with a given plate scale to a detector with given dimensions and plate scale. "
	print("Resizing image(s) to detector...")

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
		images = images[:, images.shape[height_idx]//2-detector_height_px//2:images.shape[height_idx]//2+detector_height_px//2, :]
		pad_height_top = 0
		pad_height_bottom = 0
	else:
		pad_height_top = np.floor((detector_height_px - images.shape[height_idx])/2.).astype(np.int)
		pad_height_bottom = np.ceil((detector_height_px - images.shape[height_idx])/2.).astype(np.int)

	if dest_width_px > detector_width_px:
		images = images[:, :, images.shape[width_idx]//2-detector_width_px//2:images.shape[width_idx]//2+detector_width_px//2]
		pad_width_left = 0
		pad_width_right = 0
	else: 
		pad_width_left = np.floor((detector_width_px - images.shape[width_idx])/2.).astype(np.int)
		pad_width_right = np.ceil((detector_width_px - images.shape[width_idx])/2.).astype(np.int)

	# Padding the resized images if necessary.
	images = np.pad(images, ((0, 0), (pad_height_top, pad_height_bottom), (pad_width_left, pad_width_right)), mode='constant')

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE, FIGSIZE))
		plt.subplot(1,2,1)
		plt.imshow(images_raw[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Input image')
		plt.subplot(1,2,2)
		plt.imshow(images[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Resized image')
		plt.suptitle('Resizing truth image to detector')
		plt.show()

	return np.squeeze(images)

# ##################################################################################
def getDiffractionLimitedImage(image_truth, wavelength_m, f_ratio, l_px_m, 
	detector_size_px=None,
	plotIt=False):
	" Convolve the PSF of a given telescope in a given band (J, H or K) with image_truth to simulate diffraction-limited imaging. "
	" It is assumed that the truth image has the appropriate plate scale of, but may be larger than, the detector. "
	" If the detector size is not given, then it is assumed that the input image and detector have the same dimensions. "
	print("Diffraction-limiting truth image(s)...")

	image_truth, N, height, width = getImageSize(image_truth)
	if detector_size_px != None:
		detector_height_px, detector_width_px = detector_size_px[0:2]
		if height < detector_height_px or width < detector_width_px:
			print("ERROR: the truth image must be larger than or the same size as the detector!")
			return -1		
	else:
		detector_height_px = height
		detector_width_px = width	

	# Calculating the PSF
	psf = psfKernel(wavelength_m, f_ratio, l_px_m, (2*detector_height_px, 2*detector_width_px))
	# x = np.arange(-detector_height_px//2, +detector_height_px//2 + detector_height_px%2, 1) * l_px_m
	# y = np.arange(-detector_width_px//2, +detector_width_px//2 + detector_width_px%2, 1) * l_px_m
	# X, Y = np.meshgrid(x, y)
	# # x in units of m
	# r = np.pi / wavelength_m / f_ratio * np.sqrt(np.power(X,2) + np.power(Y,2))
	# # Calculating the PSF
	# psf = np.power(2 * special.jv(1, r) / r, 2)
	# # Normalising
	# psf[psf.shape[0]/2,psf.shape[1]/2] = 1	# removing the NaN in the centre of the image
	# P = sum(psf.flatten())
	# psf /= P
	# psf = np.swapaxes(psf,0,1)
	# psf = psf.astype(np.float64)
	# Convolving the PSF and the truth image to obtain the simulated diffraction-limited image
	image_difflim = np.ndarray((N, height, width))
	for k in range(N):
		image_difflim[k] = signal.fftconvolve(psf, image_truth[k], mode='same')[detector_height_px//2 + detector_height_px%2:detector_height_px+detector_height_px//2 + detector_height_px%2 - 1, detector_width_px//2 + detector_width_px%2:detector_width_px+detector_width_px//2 + detector_width_px%2 - 1]

	if plotIt:
		plt.figure(figsize=(3*FIGSIZE, FIGSIZE))
		plt.subplot(1,3,1)
		plt.imshow(psf)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Diffraction-limited PSF of telescope')
		plt.subplot(1,3,2)
		plt.imshow(image_truth[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Truth image')
		plt.subplot(1,3,3)
		plt.imshow(image_difflim[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Diffraction-limited image')
		plt.suptitle('Diffraction-limiting image')
		plt.show()

	return np.squeeze(image_difflim)

####################################################################################################
def getSeeingLimitedImage(images, seeing_diameter_as, 
	plate_scale_as=1,
	padFactor=1,
	plotIt=False):
	"""
		 Convolve a Gaussian PSF with an input image to simulate seeing with a FWHM of seeing_diameter_as. 
	"""
	print("Seeing-limiting image(s)",end="")

	images, N, height, width = getImageSize(images)

	# Padding the source image.
	pad_ud = height // padFactor // 2
	pad_lr = width // padFactor // 2
	
	# If the image dimensions are odd, need to ad an extra row/column of zeros.
	image_padded = np.pad(images[0], ((pad_ud,pad_ud + height % 2),(pad_lr,pad_lr + width % 2)), mode='constant')
	# conv_height = image_padded.shape[0]
	# conv_width = image_padded.shape[1]
	conv_height = 2 * pad_ud + height + (height % 2)
	conv_width = 2 * pad_lr + width + (width % 2)

	# Generate a Gaussian kernel.
	kernel = np.zeros((conv_height, conv_width))
	y_as = np.arange(-conv_width//2, +conv_width//2 + conv_width%2, 1) * plate_scale_as
	x_as = np.arange(-conv_height//2, +conv_height//2 + conv_height%2, 1) * plate_scale_as
	X, Y = np.meshgrid(x_as, y_as)
	sigma = seeing_diameter_as / (2 * np.sqrt(2 * np.log(2)))
	kernel = np.exp(-(np.power(X, 2) + np.power(Y, 2)) / (2 * np.power(sigma,2)))
	kernel /= sum(kernel.flatten())
	kernel = np.pad(kernel, ((pad_ud, pad_ud + height % 2), (pad_lr, pad_lr + width % 2)), mode='constant')

	# Convolving the kernel with the image.
	image_seeing_limited = np.ndarray((N, conv_height, conv_width))
	image_seeing_limited_cropped = np.ndarray((N, height, width))

	for k in range(N):
		print('.',end="")
		image_padded = np.pad(images[k], ((pad_ud,pad_ud + height % 2),(pad_lr,pad_lr + width % 2)), mode='constant')
		image_seeing_limited[k] = signal.fftconvolve(image_padded, kernel, mode='same')
		image_seeing_limited_cropped[k] = image_seeing_limited[k,pad_ud : height + pad_ud, pad_lr : width + pad_lr]		

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE, 2*FIGSIZE))
		plt.suptitle('Seeing-limiting image')
		plt.subplot(2,2,1)
		plt.imshow(images[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Input image')
		plt.subplot(2,2,2)
		plt.imshow(kernel)
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Kernel')
		plt.subplot(2,2,3)
		plt.imshow(image_seeing_limited[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Convolved image')
		plt.subplot(2,2,4)
		plt.imshow(image_seeing_limited_cropped[0])
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Cropped, convolved image')
		plt.show()

	return np.squeeze(image_seeing_limited_cropped)

####################################################################################################
def addNoise(images,band,t_exp, 
	worstCaseSpider=False,
	plotIt=False):
	""" Add noise to an array of input images assuming an exposure time t_exp. """
	print ('Adding noise to image(s)', end="")

	# Creating an array in which to store the noisy images
	images, N, height, width = getImageSize(images)

	noisyImages = np.copy(images)

	# Getting noise parameters from the ETC.
	etc_output = exposureTimeCalc(band = band,t_exp = t_exp, worstCaseSpider = worstCaseSpider)

	# Adding noise to each image.
	for k in range(N):
		print('.', end="")
		frame_sky = np.random.poisson(lam=etc_output['N_sky'], size=(height, width))
		frame_dark = np.random.poisson(lam=etc_output['N_dark'], size=(height, width))
		frame_cryo = np.random.poisson(lam=etc_output['N_cryo'], size=(height, width))
		frame_RN = np.random.poisson(lam=etc_output['N_RN'], size=(height, width))
		noisyImages[k] += frame_sky + frame_cryo + frame_RN + frame_dark	
	print('\n')

	if plotIt:
		plt.figure(figsize=(2*FIGSIZE,FIGSIZE))
		plt.suptitle('Adding noise')
		plt.subplot(1,2,1)
		plt.imshow(images[0], vmin=min(images[0].flatten()), vmax=max(noisyImages[0].flatten()))
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Raw input image')
		plt.subplot(1,2,2)
		plt.imshow(noisyImages[0], vmin=min(images[0].flatten()), vmax=max(noisyImages[0].flatten()))
		plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
		plt.title('Noisy image')
		plt.show()

	return (np.squeeze(noisyImages), etc_output)
