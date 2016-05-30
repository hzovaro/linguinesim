############################################################################################
#
# 	File:		apd_sim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		29/05/2016
#
#	Description:
#	Simulate images of galaxies imaged using the APD detector within the ANU 2.3 m telescope.
#
#	TO DO:
#	- figure out what we're actually convolving
#		- shifting and stacking?
#		- disc-bulge decomposition? Use GALFIT for that?
#	- perform a simple convolution on two images, output the result meaningfully, 
#		shift-and-stack them.
#	- add noise to a 2D array
#	- multi-thread processing
#	- simulate tip and tilt distortion.
#
############################################################################################

import apdParameters as detector
import anu23mParameters as telescope
import cryoParameters as cryo
from apd_etc import *

import numpy as np
import numpy.fft
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.signal
import matplotlib
matplotlib.rc('image', interpolation='none', cmap = 'gray')

""" Variables """
# Imaging variables
band = 'K'
t_exp = 0.1						# exposure time
N = 10							# number of images

""" Imaging simulation """
# Add signal.
hdulist = fits.open('sample.fits')

imgs = hdulist[0].data
N = imgs.shape[0]
length = imgs.shape[1]
width = imgs.shape[2]

hdulist.close()

# Add noise. 
etc_output = exposureTimeCalc(band='H',surfaceBrightness=19,magnitudeSystem='AB',t_exp=0.1)

# The noise sigmas are expressed in units of electrons
# Need to multiply by the ADU gain to get counts
# For every pixel, get a random number drawn from a Gaussian distribution scaled up by sigma for each noise source
# and ceiling round it 

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

# img_truth = np.copy(imgs)
# img_noisy = np.copy(imgs)

# Adding noise to each image.
# for k in range(N):
# 	frame_sky = np.ones((length, width)) * etc_output['N_sky'] + np.random.randn(length, width) * etc_output['sigma_sky']
# 	frame_dark = np.ones((length, width)) * etc_output['N_dark'] + np.random.randn(length, width) * etc_output['sigma_dark']
# 	frame_cryo = np.ones((length, width)) * etc_output['N_cryo'] + np.random.randn(length, width) * etc_output['sigma_cryo']
# 	frame_RN = np.ones((length, width)) * etc_output['N_RN'] + np.random.randn(length, width) * etc_output['sigma_RN']
# 	img_noisy[k] += frame_sky + frame_cryo + frame_RN + frame_dark

# img = shiftAndStack(imgs)

# # Plotting
# plt.figure()
# plt.subplot(1,2,1)
# plt.imshow(img_truth, interpolation='nearest')
# plt.colorbar()
# plt.title('Truth image')
# plt.subplot(1,2,2)
# plt.imshow(img_noisy, interpolation='nearest')
# plt.colorbar()
# plt.title('Noisy image')
# plt.show()

# Convolutions
# Gaussian profiles
x = np.linspace(-10, +10, 60)
y = np.linspace(-10, +10, 60)
X, Y = np.meshgrid(x, y)
# I_0 = 10
# w_0 = 0.2
# kernel = I_0 * np.exp(- (np.power(X, 2) + np.power(Y, 2)) / (2 * w_0 * w_0))
psf = np.power(2 * scipy.special.jv(1,np.sqrt(np.power(X,2) + np.power(Y,2))) / (np.sqrt(np.power(X,2) + np.power(Y,2))), 2)

img_1 = np.copy(imgs[0])
img_2 = np.roll(np.copy(imgs[0]), 50, 0)
f_1 = np.fft.fft2(img_1)
f_2 = np.fft.fft2(img_2)

conv = scipy.signal.convolve2d(img_1, img_2)

# Let x and y be our two images. 
# If x is simply y shifted, then X(u,v) and Y(u,v) will be the same in frequency space, but will be shifted apart in phase space. 

plt.close('all')
plt.figure()
plt.subplot(1,2,1)
plt.imshow(img_1)
plt.title('Image 1')
plt.subplot(1,2,2)
plt.imshow(img_2)
plt.title('Image 2')

plt.figure()
plt.subplot(2,2,1)
plt.imshow(np.real(f_1))
plt.title(r'$\mathbb{Re}(\mathcal{F}(img_1))$')
plt.subplot(2,2,2)
plt.imshow(np.imag(f_1))
plt.title(r'$\mathbb{Im}(\mathcal{F}(img_1))$')
plt.subplot(2,2,3)
plt.imshow(np.real(f_2))
plt.title(r'$\mathbb{Re}(\mathcal{F}(img_2))$')
plt.subplot(2,2,4)
plt.imshow(np.imag(f_2))
plt.title(r'$\mathbb{Im}(\mathcal{F}(img_2))$')

plt.figure()
plt.imshow(np.log(conv))
plt.title('Convolution')
# plt.colorbar()
plt.show()


#########################################################################################################
# def shiftAndStack(
# 	imgs
# 	):
# 	" Shift and stack N images given in the array of N images imgs."
# 	N = imgs.shape(0)

# 	# Convolve each image w.r.t. the first.
# 	k = 1
# 	for k in range(N):
# 		# Convolve image k with image 1.
# 		# Find the peak of the convolution.

# 	# Shift-and-stack each image.

# 	return img