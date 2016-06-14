############################################################################################
#
# 	File:		apd_sim_test.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		14/06/2016
#
#	Description:
#	Testing out the functions in the apd_sim module.
#
############################################################################################

from apd_sim import *

plt.close('all')

############################################################################################

N = 5
fname = 'sample-images/tadpole.fits'
images_hst = getImages(fname)
(images, tt_idxs) = addTurbulence(images_hst[0], N, sigma=50)

image_stacked = shiftAndStack(images,plotIt=True)

plt.figure()
plt.imshow(image_stacked)
plt.colorbar()

# image_truth = np.copy(np.squeeze(images_hst[0]))

# # http://stackoverflow.com/questions/1100100/fft-based-2d-convolution-and-correlation-in-python
# width = images[0].shape[0]
# height = images[0].shape[1]
# corrs = np.ndarray((N, 2 * width - 1, 2 * height - 1))
# corr_peak_idxs = np.ndarray((N, 2))
# img_peak_idxs = np.ndarray((N, 2))

# plt.figure()
# plt.subplot(1,3,1)
# plt.imshow(images[0],origin='lower')
# plt.subplot(1,3,2)
# scat2 = plt.scatter(0.0,0.0,c='r',s=20)
# plt.subplot(1,3,3)
# scat3 = plt.scatter(0.0,0.0,c='g',s=40)

# for k in range(N):
# 	# x-correlate image k with image 1.
# 	corrs[k] = scipy.signal.fftconvolve(image_truth, images[k][::-1,::-1])
# 	corr_peak_idxs[k] = np.unravel_index(np.argmax(corrs[k]), (2 * width - 1, 2 * height - 1))
# 	img_peak_idxs[k] = np.copy(corr_peak_idxs[k])
# 	img_peak_idxs[k][0] -= width - 1
# 	img_peak_idxs[k][1] -= height - 1
# 	img_peak_idxs[k] *= -1

# 	# Plotting
# 	plt.subplot(1,3,2)
# 	plt.imshow(images[k],origin='lower')	
# 	plotcoords = np.ndarray((2))
# 	plotcoords[1] = img_peak_idxs[k,0] + width / 2
# 	plotcoords[0] = img_peak_idxs[k,1] + height / 2
# 	scat2.set_offsets(plotcoords)

# 	plt.subplot(1,3,3)
# 	plt.imshow(corrs[k],interpolation='nearest',origin='lower')
# 	corr_peak_coords = np.ndarray((2))
# 	corr_peak_coords[0] = corr_peak_idxs[k][1]
# 	corr_peak_coords[1] = corr_peak_idxs[k][0]
# 	scat3.set_offsets(corr_peak_coords)

# 	# plt.scatter([peak_idxs[k][0]], [peak_idxs[k][1]], c='r', s=20)
# 	plt.draw()
# 	plt.pause(1)


