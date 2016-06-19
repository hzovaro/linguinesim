############################################################################################
#
# 	File:		apd_sim_test.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		19/06/2016
#
#	Description:
#	Testing out the functions in the apd_sim module.
#
#	TO DO:
#	- how to estimate sigma_as for tip/tilt?
#	- get plate scales for input truth images
#	- double check the PSF stuff
#
############################################################################################

from apd_sim import *
matplotlib.rc('image', interpolation='none', cmap = 'gray')
plt.close('all')

############################################################################################

N = 5
fname = 'sample-images/tadpole.fits'
# Load the raw images.
images_hst = getImages(fname, source_plate_scale_as=0.05)
image_truth = images_hst[0]
# Diffraction-limit the image.
image_difflim = getDiffractionLimitedImage(image_truth, band='K')

# Return N images with turbulence added:
#	- randomised tip and tilt with a standard deviation of sigma_px pixels.
(images_tt, tt_idxs) = addTurbulence(image_difflim, N, sigma_px=20)
# Shift-and-stack the images.
image_stacked = shiftAndStack(images_tt, plotIt=False)

plt.figure()
plt.subplot(1,3,1)
plt.imshow(image_truth)
plt.title('Truth image')
plt.subplot(1,3,2)
plt.imshow(images_tt[0])
plt.title('Diffraction-limited, noisy image with tip-tilt')
plt.subplot(1,3,3)
plt.imshow(image_stacked)
plt.title(r'Shifted-and-stacked image, $N = %d$' % N)
# plt.colorbar()
plt.show()
