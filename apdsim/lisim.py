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
#	This file is part of lingiune-sim.
#
#	lingiune-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	lingiune-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with lingiune-sim.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################
from __future__ import division
from __future__ import print_function
from apdsim import *
import time
import pyxao

def addTipTilt(images, 
	sigma_tt_px=None,
	tt_input_idxs=None,
	N_tt=1):
	""" 
		Add turbulence to an input `truth' image. Returns N_tt copies of the input image with randomised turbulence added. 
		Just tip and tilt for now with a standard deviation of sigma_tt_px in both dimensions.

		The tip and tilt values are either random numbers drawn from a Gaussian distribution with standard deviation sigma_tt_px, or shifts specified in the input vector with shape (N, 2).
	"""

	# Tip and tilt for now	
	images, N, height, width = getImageSize(images)
	if N != 1:		
		N_tt = N # Then we add a randomised tip and tilt to each of the N input images.		
	# Otherwise, we make N_tt copies of the image and add a randomised tip and tilt to each.

	print("Adding tip/tilt to {:d} images".format(N_tt))	
	
	if not plt.is_numlike(tt_idxs):
		tt_idxs = np.ndarray((N_tt, 2))
		tt_idx = None
	else:
		sigma_tt_px = None

	# PARALLELISE
	for j in range(N_tt):
		if N == 1:			
			image = images[0]	# If the user only inputs 1 image, then we use the same image each time.
		else:
			image = images[j]

		if plt.is_numlike(tt_input_idxs):
			tt_idx = tt_input_idxs[j]			

		images_tt[j], tt_idxs[j] = addTipTilt_single(image, sigma_tt_px, tt_idx)

	return np.squeeze(images_tt), tt_idxs

####################################################################################################
def addTipTilt_single(image, 
	sigma_tt_px=None, 
	tt_idxs=None):

	if not plt.is_numlike(sigma_tt_px) and not plt.is_numlike(tt_idxs):
		print("ERROR: either sigma_tt_px OR tt_idxs must be specified!")
		raise UserWarning
	
	# Adding a randomised tip/tilt to the image
	if plt.is_numlike(sigma_tt_px):
		# If no vector of tip/tilt values is specified, then we use random numbers.
		shift_height = np.random.randn() * sigma_tt_px
		shift_width = np.random.randn() * sigma_tt_px
		tt_idxs = [shift_height, shift_width]
	else:
		# Otherwise we take them from the input vector.
		shift_height = tt_idxs[0]
		shift_width = tt_idxs[1]
	
	image_tt = shift(image, (shift_height, shift_width))

	return image_tt, tt_idxs
	

####################################################################################################
def getAoPsfs(wavelength_science_m, N_frames, psf_as_px, dt, D_out, D_in, 
	r0_500nm, v_wind_m, wind_angle_deg, elevations_m, airmass,	# Atmospheric parameters
	mode,
	N_actuators, dm_geometry, central_actuator,	# AOI: actuator in the middle
	N_lenslets, wavelength_wfs_m, wfs_geometry, central_lenslet, wfs_fratio, # fratio ~60-70 for AOI
	wave_height_px,	# Grid size in FFT (larger = better)
	plotIt = False,
	save = False, fname = "ao_psfs"
	):
	"""
		Returns a time series of PSFs (normalised by default) of a telescope with inner and outer primary mirror diameters D_in and D_out respectively in the presence of atmospheric turbulence. 
 
		The diffraction-limited PSF of the system is also returned.
	"""
	wavefrontPupil = {	
		'type':'annulus',
		'dout': D_out,
		'din' : D_in
	}

	# Wave parameters
	m_per_px = D_out / wave_height_px		# Physical mapping of wave onto primary mirror size

	# AO system parameters
	actuator_pitch_m = D_out / N_actuators
	lenslet_pitch_m = D_out / N_lenslets
	edge_radius = 1.4	
	influence_fun = 'gaussian'
	pokeStroke = 1e-7	

	# Seeing conditions
	r0_wfs = np.power((wavelength_wfs_m / 500e-9), 1.2) * r0_500nm
	r0_science = np.power((wavelength_science_m / 500e-9), 1.2) * r0_500nm

	####################################################
	# Setting up AO system
	wf_wfs = pyxao.Wavefront(wave = wavelength_wfs_m, m_per_px = m_per_px, sz = wave_height_px, pupil = wavefrontPupil)
	wf_science = pyxao.Wavefront(wave = wavelength_science_m, m_per_px = m_per_px, sz = wave_height_px, pupil = wavefrontPupil)
	wavefronts_dm = [wf_wfs, wf_science] 	# Wavefronts corrected by the DM (in a CL AO system, it's all of them!)
	wavefronts_wfs = [wf_wfs]				# Wacefronts sensed by the WFS
	psf_ix = 1		# Index in the list of wavefronts passed to the DM instance corresponding to the PSF to return

	dm = pyxao.DeformableMirror(
		wavefronts = wavefronts_dm, 
		influence_function = 'gaussian', 
		central_actuator = central_actuator, 
		actuator_pitch = actuator_pitch_m, 
		geometry = dm_geometry, 
		edge_radius = 1.4)

	wfs = pyxao.ShackHartmann(
		wavefronts = wavefronts_wfs, 
		lenslet_pitch = lenslet_pitch_m, 
		geometry = wfs_geometry, 
		central_lenslet = central_lenslet, 		
		sampling = 1)
		# fratio = wfs_fratio)

	ao = pyxao.SCFeedBackAO(dm = dm, wfs = wfs, image_ixs = psf_ix)
	
	ao.find_response_matrix()
	ao.compute_reconstructor(threshold=0.1)

	# The atmosphere is a PHASE SCREEN: it's the same at all wavelengths! We don't need to make a new atmosphere instance for each wavelength
	# If you need convincing, see minerva_red.py
	atm = pyxao.Atmosphere(sz = wave_height_px, m_per_px = m_per_px,
		elevations = elevations_m, r_0 = r0_500nm, wave_ref = 500e-9, angle_wind = wind_angle_deg,
		v_wind = v_wind_m, airmass = airmass, seed = 3)

	wf_wfs.add_atmosphere(atm)
	wf_science.add_atmosphere(atm)

	# Calculating the Nyquist oversampling factor
	psf_rad_px = np.deg2rad(psf_as_px / 3600)
	N_OS = wf_science.wave / D_out / 2 / psf_rad_px

	# Running the AO loop.
	psfs_cropped, psf_mean, psf_mean_all = ao.run_loop(dt = dt,                    # WFS photon noise.
                mode = mode,
                niter = N_frames,
                psf_ix = psf_ix,                     # Index in the list of wavefronts of the PSF you want to be returned
                plate_scale_as_px = psf_as_px,       # Plate scale of the output images/PSFs
                detector_size_px = (80,80),    # For now this is only used in plotting
                nframesbetweenplots = 5,
                plotIt = plotIt
                )

	# No turbulence.
	psf_dl = mu.centreCrop(wf_science.psf_dl(plate_scale_as_px = psf_as_px), psfs_cropped[0].shape)

	# Saving to file
	if save:
		np.savez(fname, 
			N_frames = N_frames,
			psf_atm = psf_atm,
			psf_dl = psf_dl,
			plate_scale_as_px = psf_as_px,
			N_OS = N_OS,
			dt = dt,
			wavelength_m = wavelength_m, 
			r0_500nm = r0_500nm,
			r0_science = r0_science,
			r0_lgs = r0_lgs,
			elevations = elevations,
			v_wind = v_wind,
			wind_angle = wind_angle,
			airmass = airmass,
			wave_height_px = wave_height_px,
			D_out = D_out,
			D_in = D_in
			)

	# Output the seeing-limited 
	return psfs_cropped, psf_dl

####################################################################################################
def strehl(psf, psf_dl):
	""" Calculate the Strehl ratio of an aberrated input PSF given the diffraction-limited PSF. """
	return np.amax(psf) / np.amax(psf_dl)

####################################################################################################
def shift_pp(image, img_ref_peak_idx, fsr, bid_area):
	if type(image) == list:
		image = np.array(image)	

	# Search in the bid area of the input image for the peak pixel coordinates.
	if bid_area:
		sub_image = rotateAndCrop(image_in_array = images, cropArg = bid_area)
	else:
		sub_image = image		
	img_peak_idx = np.asarray(np.unravel_index(np.argmax(sub_image), sub_image.shape))	

	# Shift the image by the relative amount.
	rel_shift_idx = (img_ref_peak_idx - img_peak_idx)
	image_shifted = shift(image, rel_shift_idx)

	peak_pixel_val = max(sub_image.flatten())	# Maximum pixel value (for now, not used)

	return image_shifted, -rel_shift_idx, peak_pixel_val

####################################################################################################
def shift_centroid(image, img_ref_peak_idx):
	if type(image) == list:
		image = np.array(image)

	img_peak_idx = _centroid(image)

	# Shift the image by the relative amount.
	rel_shift_idx = (img_ref_peak_idx - img_peak_idx)
	image_shifted = shift(image, rel_shift_idx)

	return image_shifted, -rel_shift_idx

####################################################################################################
def shift_xcorr(image, image_ref, buff, subPixelShift):
	if type(image) == list:
		image = np.array(image)

	height, width = image.shape
	corr = signal.fftconvolve(image_ref, image[::-1,::-1], 'same')
	corr /= max(corr.flatten())	# The fitting here does not work if the pixels have large values!
	
	if subPixelShift: 
		# Fitting a Gaussian.
		Y, X = np.mgrid[-(height-2*buff)/2:(height-2*buff)/2, -(width-2*buff)/2:(width-2*buff)/2]		
		p_init = models.Gaussian2D(x_stddev=1.,y_stddev=1.)
		fit_p = fitting.LevMarLSQFitter()
		p_fit = fit_p(p_init, X, Y, corr[buff:height-buff, buff:width-buff])		
		rel_shift_idx = (p_fit.y_mean.value, p_fit.x_mean.value)	# NOTE: the indices have to be swapped around here for some reason!		
	else:
		rel_shift_idx = np.unravel_index(np.argmax(corr), corr.shape)
		rel_shift_idx = (rel_shift_idx[0] - height/2, rel_shift_idx[1] - width/2)
	
	image_shifted = shift(image, rel_shift_idx)	

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
	print("Applying Lucky Imaging technique '{}' to input series of {:d} images...".format(li_method, N))
	
	# For each of these functions, the output must be of the form 
	#	image_shifted, rel_shift_idxs	
	if li_method == 'xcorr':
		shift_fun = partial(shift_xcorr, image_ref=image_ref, buff=buff, subPixelShift=subPixelShift)
	
	elif li_method == 'peak_pixel':
		# Determining the reference coordinates.
		if bid_area:			
			sub_image_ref = rotateAndCrop(image_in_array = image_ref, cropArg = bid_area)
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
		pool = ThreadPool()
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
		print("Elapsed time for {:d} {}-by-{} images in {} mode: {:.5f}".format(N, image_ref.shape[0], image_ref.shape[1], mode, (toc-tic)))

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
	verbose=True):
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
	if image_ref != None and type(image_ref.flatten()[0]) != np.float64:
		image_ref = image_ref.astype(np.float64)

	# Checking image dimensions.
	if len(images.shape) > 4:
		print("WARNING: for now, please only input a 3D array of images to shift and stack! I'm only going to operate on the first set of images...")
		images = np.squeeze(images[0])	
	
	if len(images.shape) == 3:
		if N and N > images.shape[0]:
			print("ERROR: if specified, N must be equal to or less than the length of the first dimension of the images array.")
			raise UserWarning
		if image_ref == None:
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