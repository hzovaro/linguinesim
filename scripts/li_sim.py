############################################################################################
#
# 	File:		li_sim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		29/06/2016
#
#	Description:
#	Generating some simulated images of the lucky imaging process. Not worrying about 
#	resizing to the detector for now.
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
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
############################################################################################
from __future__ import division
from apd_sim import *
plt.close('all')

N_exp = 5
t_exp = 10e-3
band = 'K'
sigma_tt_px = 20

" ARP 194 "
# fname = 'sample-images/arp194/arp194.fits'
# angle = 55+90.
# cropIdx = (1000,500,2500,1500)

" ARP 274 "
fname = 'sample-images/arp274/hlsp_heritage_hst_wfpc2_arp274_f814w_v1_sci.fits'
angle = 10.
cropIdx = (500,300,2700,1500)

wfpc2_plate_scale_as_1 = 1.38889E-05 * 3600.
wfpc2_plate_scale_as_2 = 2.77778E-05 * 3600.
plate_scale_as = wfpc2_plate_scale_as_1
seeing_diameter_as = 1.5

# 1. Load the truth image.
fitsFileData = getRawImages(fname)
images_hst = fitsFileData[0]

# 2. Resize and crop to make it look pretty.
images_hst = rotateAndCrop(images_hst, angle=angle, cropIdx=cropIdx)
image_diffraction_limited = np.fliplr(images_hst)	# Flipping because that's what Rob's one looks like

detector_height_px = image_diffraction_limited.shape[0]
detector_width_px = image_diffraction_limited.shape[1]

# TODO: resizing the image does something weird to the pixel values
# image_diffraction_limited_small = resizeImagesToDetector(images_hst, 
# 	source_plate_scale_as = plate_scale_as, 
# 	dest_plate_scale_as = plate_scale_as * 2, 
# 	dest_detector_size_px = (detector_height_px/2, detector_width_px/2),
# 	plotIt=True)

# 3. Seeing-limit the image by convolving with a Gaussian with a specified FWHM.
# image_seeing_limited = getSeeingLimitedImage(image_diffraction_limited_small, seeing_diameter_as, plate_scale_as * 2, plotIt=True)
image_seeing_limited = getSeeingLimitedImage(image_diffraction_limited, seeing_diameter_as, plate_scale_as, plotIt=True)

# 4. Generate N images with randomised tip and tilt.
images_tt, tt_idx = addTurbulence(image_seeing_limited, 
	N = N_exp, 
	sigma_tt_px = sigma_tt_px)

# 5. Shift-and-stack the images.
image_stacked = shiftAndStack(images_tt, image_ref=image_seeing_limited, plotIt=False)

# Plotting
subplot_x = 3
subplot_y = 1
figsize = 5.
plt.figure(figsize=(subplot_x * figsize, subplot_y * figsize))
plt.subplot(subplot_y,subplot_x,1)
plt.imshow(image_diffraction_limited, norm=LogNorm())
plt.title('Diffraction limited image')
plt.subplot(subplot_y,subplot_x,2)
# plt.imshow(np.log(image_diffraction_limited_small))
plt.imshow(np.log(image_seeing_limited))
plt.title('Seeing limited image, seeing = %.1f"' % seeing_diameter_as)
plt.subplot(subplot_y,subplot_x,3)
plt.imshow(image_stacked, norm=LogNorm())
plt.title('Shifted and stacked image')
plt.show()

plt.figure(figsize=(7.5,7.5))
plt.imshow(np.log(image_diffraction_limited))
plt.title('HST diffraction limited image')
plt.figure(figsize=(7.5,7.5))
plt.imshow(image_seeing_limited, norm=LogNorm())
plt.title('Seeing limited image, seeing = %.1f"' % seeing_diameter_as)
plt.show()