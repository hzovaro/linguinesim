#########################################################################################################
#
# 	File:		apd_sim_test.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Testing out the functions in the apdsim package.
#
#	Copyright (C) 2016 Anna Zovaro
#
#########################################################################################################
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
#########################################################################################################
#
#	TO DO:
#	- in-quadrature PSF: figure out sigma for 
#	- simulate realistic magnitudes for target objects
#	- double check the PSF stuff
#
#########################################################################################################
from __future__ import division
from apdsim import *
plt.close('all')

N_exp = 200
t_exp = 20e-3
band = 'K'

targname = "ARP 194"
fname = 'sample-images/arp194/arp194.fits'
# Cropping & rotating the raw input image to make it more manageable (and prettier)
angle = 55+90.
# TODO: 'centre' the image here, before we pass it any further. So define a source image width and height (in px)
# and then define the crop indices relative to the centre coordinates.
crop_left = 1100
crop_right = 2300
crop_top = 400
crop_bottom = 1400
cropIdx = (crop_left,crop_top,crop_right,crop_bottom)

# targname = "ARP 274"
# fname = 'sample-images/arp274/hlsp_heritage_hst_wfpc2_arp274_f814w_v1_sci.fits'
# # Cropping & rotating the raw input image to make it more manageable (and prettier)
# angle = 10.
# crop_left = 500
# crop_right = 2700
# crop_top = 300
# crop_bottom = 1500
# cropIdx = (crop_left,crop_top,crop_right,crop_bottom)

# HST & system plate scales
wfpc2_plate_scale_as_1 = 1.38889E-05 * 3600
wfpc2_plate_scale_as_2 = 2.77778E-05 * 3600
sys_plate_scale_as = telescope.plate_scale_as_m * detector.l_px_m
f_ratio = telescope.f_ratio

# Tip/tilt parameters
pad_tt = 50			# Pad the input image before applying tip/tilt so that we don't get zeros at the edges after tip and tilt are applied.
seeing_as = 2		# Seeing disc diameter (FWHM) in arcsec
sigma_tt_px = seeing_as / 2 / sys_plate_scale_as

# Get the raw images.
fitsFileData = getRawImages(fname)
images_hst = fitsFileData[0]
images_hst, N = getImageSize(images_hst)[0:2]
# if len(images_hst.shape) > 2:
# 	N = images_hst.shape[0]
# else:
# 	N = 1

# Rotate and crop.
images_hst = rotateAndCrop(images_hst, angle=angle, cropArg=cropIdx, plotIt=True)
# images_hst += np.abs(min(images_hst.flatten()))
# images_hst = np.fliplr(images_hst)	# Flipping because that's what Rob's one looks like

# Resize to the detector.
images_resized = resizeImagesToDetector(images_hst, 
	source_plate_scale_as = wfpc2_plate_scale_as_1, 
	dest_plate_scale_as = sys_plate_scale_as, 
	dest_detector_size_px = (detector.height_px + 2 * pad_tt, detector.width_px + 2 * pad_tt),
	plotIt=True)

############################################################################################
# EXPERIMENTAL: adjusting for surface brightness
# NOTE: do this after resizing
# 1. normalise the image
if N > 1:
	for k in range(N):
		maxval = max(images_resized[k].flatten())
		images_resized[k] /= maxval
else:
	maxval = max(images_resized.flatten())
	images_resized /= maxval

# 2. Multiply by the flux * t_exp
flux = getPhotonFlux(surfaceBrightness = 12, 
		wavelength_eff = telescope.filter_bands_m[band][0], 
		bandwidth = telescope.filter_bands_m[band][1], 
		plate_scale_as = sys_plate_scale_as, 
		A_tel = telescope.A_collecting, 
		tau = telescope.tau,
		qe = detector.qe,
		gain_av = detector.gain_av,
		magnitudeSystem = 'AB')

images_resized *= t_exp * flux

############################################################################################

# Diffraction-limit the truth (HST) image.
image_difflim_padded = getDiffractionLimitedImage(images_resized, 
	f_ratio = f_ratio,
	detector_size_px = (detector.height_px + 2 * pad_tt, detector.width_px + 2 * pad_tt),
	l_px_m = detector.l_px_m,
	wavelength=telescope.filter_bands_m[band][0],
	plotIt=False)

image_difflim = rotateAndCrop(image_difflim_padded, angle=0., 
	cropArg = pad_tt)

# Return N images with tip/tilt added.
images_tt, tt_idxs = addTurbulence(image_difflim_padded, 
	N_tt = N_exp, 
	crop_tt = pad_tt,
	sigma_tt_px = sigma_tt_px)

# Add noise.
images_noisy, etc_output = addNoise(images_tt, 
	band = band, 
	t_exp = t_exp)
# images_noisy = np.copy(images_tt)

# Shift-and-stack the images.
N_vals = np.array([50])
images_stacked = np.ndarray((len(N_vals), images_tt.shape[1], images_tt.shape[2]))
for k in range(len(N_vals)):
	images_stacked[k] = shiftAndStack(images_noisy, image_ref=image_difflim, N=N_vals[k], plotIt=False)

# Plotting.
# vmin = min(images_hst.flatten())
plt.figure()
if N == 1:
	plt.imshow(images_hst, norm=LogNorm())
else:
	plt.imshow(images_hst[0], norm=LogNorm())
plt.title(targname + ' raw HST image')
plt.colorbar(fraction=0.046, pad=0.04)

plt.figure()
plt.imshow(images_noisy[0], norm=LogNorm())
# plt.imshow(images_noisy[0], norm=LogNorm(vmin = min(images_noisy.flatten()), vmax = max(image_stacked.flatten())))
# plt.imshow(images_noisy[0], vmin = min(images_noisy.flatten()), vmax = max(image_stacked.flatten()))
plt.title('Diffraction-limited noisy image with tip-tilt')
plt.colorbar(fraction=0.046, pad=0.04)

vmax = max(images_stacked[-1].flatten())
vmin = min(images_stacked[0].flatten())
for k in range(len(N_vals)):
	plt.figure()
	# plt.imshow(images_stacked[k], norm=LogNorm(vmin=vmin,vmax=vmax))
	plt.imshow(images_stacked[k], vmin=vmin,vmax=vmax)
	plt.title(r'Shifted-and-stacked image, $N = %d$' % N_vals[k])
	plt.colorbar(fraction=0.046, pad=0.04)
plt.show()
