############################################################################################
#
# 	File:		psf_in_quadrature.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Demonstrating that PSFs and Gaussian blur can be 'applied' successively in quadrature.
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
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
############################################################################################
from __future__ import division
from apdsim import *
from cosmo_calc import *
plt.close('all')

imsize = 501
image_truth = np.zeros((imsize, imsize))
image_truth[imsize//2+1, imsize//2+1] = 1

sigma_1 = 25
sigma_2 = 30
sigma_3 = np.sqrt(sigma_1*sigma_1 + sigma_2*sigma_2)

image_1 = getSeeingLimitedImage(images = image_truth, seeing_diameter_as = sigma_1, plotIt = False)
image_2 = getSeeingLimitedImage(images = image_truth, seeing_diameter_as = sigma_2, plotIt = False)

image_3 = getSeeingLimitedImage(images = image_1, seeing_diameter_as = sigma_2, plotIt = False)

image_4 = getSeeingLimitedImage(images = image_truth, seeing_diameter_as = sigma_3, plotIt = False)

image_1 /= max(image_1.flatten())
image_2 /= max(image_2.flatten())
image_3 /= max(image_3.flatten())
image_4 /= max(image_4.flatten())

# plt.figure()
# plt.plot(image_1[imsize//2,:], 'm')
# plt.plot(image_2[imsize//2,:], 'k')
# plt.plot(image_4[imsize//2,:], 'b')
# plt.plot(image_3[imsize//2,:], 'r')
# plt.show()

# Same but in the diffraction limit...
efl = 20
wavelength_1 = 500e-9
wavelength_2 = 800e-9
D_1 = 0.1
D_2 = 0.05
fwhm_1 = wavelength_1 / D_1
fwhm_2 = wavelength_2 / D_2
f_ratio_1 = efl / D_1
f_ratio_2 = efl / D_2

fwhm_3 = np.sqrt(fwhm_2*fwhm_2 - fwhm_1*fwhm_1)
wavelength_3 = wavelength_1
D_3 = wavelength_3 / fwhm_3
f_ratio_3 = efl / D_3

image_1 = getDiffractionLimitedImage(image_truth, l_px_m = 5e-6, f_ratio = f_ratio_1, wavelength_m = wavelength_1, detector_size_px = (imsize, imsize),plotIt = False)

image_2 = getDiffractionLimitedImage(image_truth, f_ratio = f_ratio_2, detector_size_px = (imsize, imsize), l_px_m = 5e-6, wavelength_m = wavelength_2, plotIt = False)

# image_3 = getDiffractionLimitedImage(image_1, f_ratio = f_ratio_3, detector_size_px = (imsize, imsize), l_px_m = 5e-6, wavelength_m = wavelength_3, plotIt = False)

image_4 = getDiffractionLimitedImage(
	image_1, 
	f_ratio = f_ratio_2, 
	f_ratio_in = f_ratio_1, 
	detector_size_px = (imsize, imsize), 
	l_px_m = 5e-6, 
	wavelength_m = wavelength_2, 
	wavelength_in_m = wavelength_1, 
	plotIt = True)
	
image_1 /= max(image_1.flatten())
image_2 /= max(image_2.flatten())
image_3 /= max(image_3.flatten())
image_4 /= max(image_4.flatten())

plt.figure()
plt.plot(image_1[imsize//2,:], 'm', label='Input image')
plt.plot(image_2[imsize//2,:], 'k', label='Ideal output image')
plt.plot(image_4[imsize//2,:], 'r', label='Convolved output image')
plt.xlabel('x (pixels)')
plt.xlim([0, imsize-1])
plt.legend(loc='southwest')
plt.show()