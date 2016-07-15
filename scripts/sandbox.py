############################################################################################
#
# 	File:		sandbox.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A 'sandbox' for testing various features. Play safe kids!
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
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
############################################################################################
from __future__ import division
from apdsim import *
from astropy.modeling import models, fitting
plt.close('all')

height = 256
width = 320
I_0 = 1e3
sigma = 10
x_0 = 10.1
y_0 = 30.1

# Sigma?
# Background noise levels?

Y, X = np.mgrid[-height/2:height/2, -width/2:width/2]
# Y, X = np.mgrid[:height, :width]

# Z = I_0 * np.exp(-((X - x_0)**2 + (Y - y_0)**2) / (2 * sigma**2) )
image_1 = getStarField(m = 10, coords = np.array((height/2, width/2)), t_exp = 10e-3, A_tel = telescope.A_collecting, f_ratio = telescope.f_ratio*20, l_px_m = detector.l_px_m, detector_size_px = detector.size_px, magnitudeSystem = 'AB', band = 'K', plotIt = False, tau = telescope.tau, gain = detector.gain, qe = detector.qe, detectorSaturation = detector.saturation)[1]
image_2 = shift(image_1, (x_0, y_0))
# image_2 = np.copy(image_1)
corr = signal.fftconvolve(image_1, image_2[::-1,::-1],'same')
corr = corr / max(corr.flatten())
# corr = corr[height/2-1:3*height/2-1,width/2-1:3*width/2-1]

p_init = models.Gaussian2D(x_stddev=1.,y_stddev=1.)
fit_p = fitting.LevMarLSQFitter()
p_fit = fit_p(p_init, X, Y, corr)

max_idx = np.array((p_fit.x_mean.value, p_fit.y_mean.value))
print("Indices of maximum in Gaussian fit relative to image centre: (%5.2f,%5.2f)" % (max_idx[0], max_idx[1]))

plt.figure(figsize=(2*FIGSIZE,FIGSIZE))
plt.subplot(1,2,1)
plt.imshow(corr)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title('Cross-correlation')
plt.subplot(1,2,2)
plt.imshow(p_fit(X,Y))
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title('Model fit')
plt.show()