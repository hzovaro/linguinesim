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

def addTurbulence(images, 
	sigma_tt_px=None,
	tt_idxs=None,
	N_tt=1,
	crop_tt=None):
	""" 
		Add turbulence to an input `truth' image. Returns N_tt copies of the input image with randomised turbulence added. 
		Just tip and tilt for now with a standard deviation of sigma_tt_px in both dimensions.

		The tip and tilt values are either random numbers drawn from a Gaussian distribution with standard deviation sigma_tt_px, or shifts specified in the input vector with shape (N, 2).
	"""

	# Tip and tilt for now	
	images, N, height, width = getImageSize(images)
	print("Adding tip/tilt to {:d} images".format(N))

	if not plt.is_numlike(sigma_tt_px) and not plt.is_numlike(tt_idxs):
		print("ERROR: either sigma_tt_px OR tt_idxs must be specified!")
		raise UserWarning

	if N != 1:
		# Then we add a randomised tip and tilt to each of the N input images.
		N_tt = N
	# Otherwise, we make N_tt copies of the image and add a randomised tip and tilt to each.

	# Output array of images
	if not plt.is_numlike(crop_tt):
		images_tt = np.ndarray((N_tt, height, width))
	else:
		if len(crop_tt.shape) == 0:
			images_tt = np.ndarray((N_tt, height - 2 * crop_tt, width - 2 * crop_tt))	
		else:
			images_tt = np.ndarray((N_tt, height - 2 * crop_tt[0], width - 2 * crop_tt[1]))

	# Array to hold the tip/tilt offsets
	if not plt.is_numlike(tt_idxs):
		tt_idxs = np.ndarray((N_tt, 2))
	
	# Adding a randomised tip/tilt to each of N_tt images
	for j in range(N_tt):
		if N == 1:
			image = images[0]
		else:
			image = images[j]

		if plt.is_numlike(sigma_tt_px):
			# If no vector of tip/tilt values is specified, then we use random numbers.
			shift_height = np.random.randn() * sigma_tt_px
			shift_width = np.random.randn() * sigma_tt_px
			tt_idxs[j] = [shift_height, shift_width]
		else:
			# Otherwise we take them from the input vector.
			shift_height = tt_idxs[j,0]
			shift_width = tt_idxs[j,1]
		
		image_tt = shift(image, (shift_height, shift_width))

		# Cropping the image if necessary
		if crop_tt == None:
			images_tt[j] = image_tt
		else:
			images_tt[j] = rotateAndCrop(image_tt, angle=0., cropArg=crop_tt)
	return np.squeeze(images_tt), np.squeeze(tt_idxs)

################################################################################################################
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

##########################################################################################
def strehl(psf, psf_dl):
	""" Calculate the Strehl ratio of an aberrated input PSF given the diffraction-limited PSF. """
	return np.amax(psf) / np.amax(psf_dl)

####################################################################################################
def xcorrShiftAndStack(images, 
	image_ref=None,
	N=None, 
	subPixelShift=None,			# Use a 2D Gaussian fit to locate the maximm in the cross-correlation output.
	buff=None,					# Used in the 2D Gaussian fit.
	showAnimatedPlots=False,	# Show an animated plot window of the shifting-and-stacking process.
	plotIt=False):
	""" 
		Shift and stack the images given in the 3-dimensional array of N images. 
	"""
	print("Applying Lucky Imaging technique: cross-correlation shift-and-stack method")

	# Input argument validation
	if not buff:
		buff = 25
	if subPixelShift == None:
		subPixelShift = True

	images, image_ref, N = _li_error_check(images, image_ref, N)	
	height = images[0].shape[0]
	width = images[0].shape[1]
	image_stacked = np.copy(image_ref)		# array to hold shifted-and-stacked image	
	
	corrs = np.ndarray((N, height, width))	# Array to hold x-correlation results
	corr_peak_idxs = np.ndarray((N, 2))		# indices of the peak value in the x-correlation
	rel_shift_idxs = np.ndarray((N, 2))		# shift in x and y computed from the x-correlation

	for k in range(N):
		# Cross-correlate image k with the reference image to find the tip and tilt.
		corrs[k] = signal.fftconvolve(image_ref, images[k][::-1,::-1], 'same')
		corrs[k] /= max(corrs[k].flatten())	# The fitting here does not work if the pixels have large values!
		if subPixelShift: 
			# Fitting a Gaussian.
			# buff = 25
			Y, X = np.mgrid[-(height-2*buff)/2:(height-2*buff)/2, -(width-2*buff)/2:(width-2*buff)/2]		
			p_init = models.Gaussian2D(x_stddev=1.,y_stddev=1.)
			fit_p = fitting.LevMarLSQFitter()
			p_fit = fit_p(p_init, X, Y, corrs[k][buff:height-buff, buff:width-buff])
			# NOTE: the indices have to be swapped around here for some reason!
			corr_peak_idxs[k] = np.array((p_fit.y_mean.value, p_fit.x_mean.value))
			rel_shift_idxs[k][0] = corr_peak_idxs[k][0]
			rel_shift_idxs[k][1] = corr_peak_idxs[k][1]

			# if plotIt:
			# 	fig = mu.newfigure(1,2)
			# 	plt.subplot(1,2,1)
			# 	plt.imshow(corrs[k])
			# 	mu.colourbar()
			# 	plt.title('Cross-correlation')
			# 	plt.subplot(1,2,2)
			# 	plt.imshow(p_fit(X,Y))
			# 	mu.colourbar()
			# 	plt.title('Model fit')
			# 	plt.show()	
			# 	# pdb.set_trace()
			# 	plt.close(fig)

		else:
			corr_peak_idxs[k] = np.unravel_index(np.argmax(corrs[k]), corrs[k].shape)
			rel_shift_idxs[k][0] = corr_peak_idxs[k][0]
			rel_shift_idxs[k][1] = corr_peak_idxs[k][1]

		# Shift-and-stack the images.
		image_stacked += shift(images[k], (rel_shift_idxs[k][0], rel_shift_idxs[k][1]))

		# Plotting
		if showAnimatedPlots:
			if k == 0:
				mu.newfigure(1,3)
				plt.suptitle('Cross-correlation Lucky Imaging output')
				plt.subplot(1,3,1)
				plt.title('Mean-combined shifted-and-stacked image')
				plt.subplot(1,3,2)
				plt.title('Single exposure')
				scat2 = plt.scatter(0.0,0.0,c='r',s=20)
				plt.subplot(1,3,3)
				plt.title('Cross-correlation')
				scat3 = plt.scatter(0.0,0.0,c='g',s=40)

			plt.subplot(1,3,1)
			plt.imshow(image_stacked)

			plt.subplot(1,3,2)
			plt.imshow(images[k])	
			plotcoords = np.ndarray((2))
			plotcoords[1] = rel_shift_idxs[k,0] + width / 2
			plotcoords[0] = rel_shift_idxs[k,1] + height / 2
			scat2.set_offsets(plotcoords)

			plt.subplot(1,3,3)
			plt.imshow(corrs[k])
			corr_peak_coords = np.ndarray((2))
			corr_peak_coords[0] = corr_peak_idxs[k][1]
			corr_peak_coords[1] = corr_peak_idxs[k][0]
			scat3.set_offsets(corr_peak_coords)

			# plt.scatter([peak_idxs[k][0]], [peak_idxs[k][1]], c='r', s=20)
			plt.draw()
			plt.pause(0.5)
	# Average-combining the images 
	image_stacked /= N + 1	# Need to add 1 because of the reference image

	if plotIt:
		mu.newfigure(1,2)
		plt.suptitle('Cross-correlation Lucky Imaging output')
		plt.subplot(1,2,1)
		# plt.imshow(images[0], vmin=min(images[0].flatten()), vmax=max(image_stacked.flatten()))
		plt.imshow(images[0])
		plt.title('Single exposure')
		mu.colourbar()
		plt.subplot(1,2,2)
		# plt.imshow(image_stacked, vmin=min(images[0].flatten()), vmax=max(image_stacked.flatten()))
		plt.imshow(image_stacked)
		plt.title('Mean-combined shifted-and-stacked image')
		mu.colourbar()
		plt.show()

	
	return image_stacked, -rel_shift_idxs

####################################################################################################
def fourierShiftAndStack(images, 
	image_ref = None, 
	plotIt = False):
	"""
		A Fourier-space method of Lucky Imaging.
		Credit: Garrel et al. (2011)
	"""
	print("Applying Lucky Imaging technique: Fourier-space method")
	images, image_ref, N = _li_error_check(images, N)	
	height = images[0].shape[0]
	width = images[0].shape[1]

	return image_stacked, other_data

####################################################################################################
def peakPixelShiftAndStack(images, 
	N = None,			# How many images we want to iterate through
	image_ref = None, 	# The reference image. By default the first image in the input images array
	bid_area = None,	# Subwindow in which to calculate the offsets. By default the entire image.
	fsr = 1,			# Frame selection rate.
	plotIt = False,
	showAnimatedPlots = False):	
	"""
		A Lucky Imaging method which uses the brightest pixel in each image to align them.
		Credit: Aspin et al. (1997)
	"""
	print("Applying Lucky Imaging technique: peak pixel shift-and-stack method, FSR = {:.2f}%".format(fsr*100))
	images, image_ref, N = _li_error_check(images, image_ref, N)	
	height = images[0].shape[0]
	width = images[0].shape[1]
	image_stacked = np.copy(image_ref)

	# 1. Search in the bid area of the reference image for the peak pixel coordinates.
	if not bid_area:
		# If no bid area is specified, then we just search the whole image.
		sub_image_ref = np.copy(image_ref)
		sub_images = np.copy(images)
	else:
		sub_image_ref = rotateAndCrop(image_in_array = image_ref, cropArg = bid_area)	# (left, upper, right, lower)
		sub_images = rotateAndCrop(image_in_array = images, cropArg = bid_area)	# (left, upper, right, lower)

	# 2. Iterate through each image and repeat.
	img_ref_peak_idx = np.asarray(np.unravel_index(np.argmax(sub_image_ref), sub_image_ref.shape)) # coordinates to align all the images to.
	img_peak_idxs = np.zeros((N, 2)) 	# Coordinates of the peak pixel in the subwindow.
	rel_shift_idxs = np.zeros((N, 2))	# Coordinates by which to shift the image (relative to image_ref)

	# Finding the peak value in each image.
	# pdb.set_trace()
	img_peak_vals = np.amax(sub_images,(1,2))
	sorted_idx = np.argsort(img_peak_vals)[::-1]
	final_idx = np.ceil(N * fsr)
	N_frames_used = np.ceil(N * fsr) + 1	# The +1 includes the reference image!

	firstLoop = True
	for k in sorted_idx[:final_idx]:
		img_peak_idxs[k] = np.asarray(np.unravel_index(np.argmax(sub_images[k]), sub_images[k].shape))

		# 3. Shift the image by the relative amount.
		rel_shift_idxs[k] = (img_ref_peak_idx - img_peak_idxs[k])
		image_stacked += shift(images[k], (rel_shift_idxs[k][0], rel_shift_idxs[k][1]))

		# Plotting
		if showAnimatedPlots:
			if firstLoop:
				mu.newfigure(1,2)
				plt.suptitle('Peak pixel Lucky Imaging output')
				plt.subplot(1,2,2)
				plt.title('Mean-combined shifted-and-stacked image')
				plt.subplot(1,2,1)
				plt.title('Single exposure')
				scat2 = plt.scatter(0.0,0.0,c='r',s=20)
				scat3 = plt.scatter(0.0,0.0,c='g',s=20)
				firstLoop = False

			plt.subplot(1,2,2)
			plt.imshow(image_stacked)

			plt.subplot(1,2,1)
			plt.imshow(images[k])	
			plotcoords = np.ndarray((2))
			plotcoords[1] = rel_shift_idxs[k,0] + width / 2
			plotcoords[0] = rel_shift_idxs[k,1] + height / 2
			scat2.set_offsets(plotcoords)
			scat3.set_offsets(np.asarray((img_peak_idxs[k,1], img_peak_idxs[k,0])))

			plt.draw()
			plt.pause(1)
	# Mean combining
	image_stacked /= N_frames_used

	if plotIt:
		mu.newfigure(1,2)
		plt.suptitle('Peak pixel Lucky Imaging output, FSR = {:.2f}%, {:d} frames used'.format(fsr * 100, int(N_frames_used)))
		plt.subplot(1,2,1)
		plt.imshow(sub_images[0])
		plt.title('Single exposure')
		plt.subplot(1,2,2)
		plt.imshow(image_stacked)
		plt.title('Mean-combined shifted-and-stacked image')
		plt.show()

	return image_stacked, -rel_shift_idxs

####################################################################################################
def centroidShiftAndStack(images, 
	N = None,			# How many images we want to iterate through
	image_ref = None, 	# The reference image. By default the first image in the input images array
	bid_area = None,	# Subwindow in which to calculate the offsets.
	plotIt = False,
	showAnimatedPlots = False):	
	"""
		A Lucky Imaging method which uses the brightest pixel in each image to align them.
		Credit: Aspin et al. (1997)
	"""
	print("Applying Lucky Imaging technique: centroiding shift-and-stack method")
	images, image_ref, N = _li_error_check(images, image_ref, N)	
	height = images[0].shape[0]
	width = images[0].shape[1]
	image_stacked = np.copy(image_ref)

	# 1. Search in the bid area of the reference image for the peak pixel coordinates.
	if not bid_area:
		# If no bid area is specified, then we just search the whole image.
		sub_image_ref = np.copy(image_ref)
		sub_images = np.copy(images)
	else:
		sub_image_ref = rotateAndCrop(image_in_array = image_ref, cropArg = bid_area)	# (left, upper, right, lower)
		sub_images = rotateAndCrop(image_in_array = images, cropArg = bid_area)	# (left, upper, right, lower)

	# 2. Iterate through each image and repeat.
	img_ref_peak_idx = _centroid(image_ref)
	img_peak_idxs = np.zeros((N, 2)) 	# Coordinates of the peak pixel in the subwindow.
	rel_shift_idxs = np.zeros((N, 2))	# Coordinates by which to shift the image (relative to image_ref)
	for k in range(N):
		img_peak_idxs[k] = _centroid(images[k])
		# 3. Shift the image by the relative amount.
		rel_shift_idxs[k] = (img_ref_peak_idx - img_peak_idxs[k])
		image_stacked += shift(images[k], (rel_shift_idxs[k][0], rel_shift_idxs[k][1]))

		# Plotting
		if showAnimatedPlots:
			if k == 0:
				mu.newfigure(1,2)
				plt.suptitle('Centroiding Lucky Imaging output')
				plt.subplot(1,2,2)
				plt.title('Mean-combined shifted-and-stacked image')
				plt.subplot(1,2,1)
				plt.title('Single exposure')
				scat2 = plt.scatter(0.0,0.0,c='r',s=20)
				scat3 = plt.scatter(0.0,0.0,c='g',s=20)

			plt.subplot(1,2,2)
			plt.imshow(image_stacked)

			plt.subplot(1,2,1)
			plt.imshow(images[k])	
			plotcoords = np.ndarray((2))
			plotcoords[1] = rel_shift_idxs[k,0] + height / 2
			plotcoords[0] = rel_shift_idxs[k,1] + width / 2
			scat2.set_offsets(plotcoords)
			scat3.set_offsets(np.asarray((img_peak_idxs[k,1], img_peak_idxs[k,0])))

			plt.draw()
			# pdb.set_trace()
			plt.pause(0.5)
	# Mean combining
	image_stacked /= N + 1	# Need to add 1 because of the reference image

	if plotIt:
		mu.newfigure(1,2)
		plt.suptitle('Centroiding Lucky Imaging output')
		plt.subplot(1,2,1)
		plt.imshow(sub_images[0])
		plt.title('Single exposure')
		plt.subplot(1,2,2)
		plt.imshow(image_stacked)
		plt.title('Mean-combined shifted-and-stacked image')
		plt.show()

	return image_stacked, -rel_shift_idxs


####################################################################################################
def luckyImaging(images, type, 
	fsr = 1,			# for peak pixel method
	bid_area = None,	# for peak pixel method
	N = None,
	image_ref = None,
	plotIt = False,
	showAnimatedPlots = False):
	""" 
		Apply a Lucky Imaging (LI) technique to a sequence of images stored in the input array images. 
		The type of LI technique used is specified by input string type and any additional arguments which may be required are given in vararg.
	"""
	if type == 'xcorr':
		image_lucky, shift_idxs = xcorrShiftAndStack(images = images, image_ref = image_ref, plotIt=plotIt, showAnimatedPlots=showAnimatedPlots)
	elif type == 'fourier':
		image_lucky = fourierLuckyImaging(images = images, image_ref = image_ref, plotIt=plotIt)
	elif type == 'peak_pixel':
		image_lucky, shift_idxs = peakPixelShiftAndStack(images = images, image_ref = image_ref, fsr = fsr, bid_area = bid_area, plotIt=plotIt, showAnimatedPlots=showAnimatedPlots)
	elif type == 'centroid':
		image_lucky, shift_idxs = centroidShiftAndStack(images = images, image_ref = image_ref, plotIt=plotIt, showAnimatedPlots=showAnimatedPlots)
	return image_lucky, shift_idxs

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