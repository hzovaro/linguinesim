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
from cosmo_calc import *
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
plt.close('all')
 
# First: 2D case. Then move onto 3D case.

# PROBLEM: the Airy function is given in units of intensity I(x, y) (in photons/s/m^2). We, however, need the photon count rate (in photons/s). So we need to integrate the Airy function over each pixel to get the count rate. 
# This needs to be done in 2D.

# Using a lambda function for the Airy disc:

# 1. Make a function for the integrand. 
wavelength_m = FILTER_BANDS_M['K'][0]
f_ratio = telescope.f_ratio
# Total intensity (W)
P_0 = 1
# Central intensity (W m^-2)
I_0 = P_0 * np.pi / 4 / wavelength_m / wavelength_m / f_ratio / f_ratio


# Calculating the PSF
# x, y will be in units of metres
r = lambda x, y: np.pi / wavelength_m / f_ratio * np.sqrt(np.power(x,2) + np.power(y,2))
I = lambda x, y : np.power((2 * special.jv(1, r(x,y)) / r(x,y)), 2) * I_0 

# 2. Calculate the sum of the integral for each pixel.
# Sanity check: over the infinite plane. Should get back P_0.
# P = integrate.nquad(I, [[-100*l_px_m, +100*l_px_m], [-100*l_px_m, +100*l_px_m]])
# print P

# Evaluating I for a grid.
l_px_m = detector.l_px_m
oversampleFactor = 4
sz_x = 3
sz_y = 3
detector_height_px = sz_x * oversampleFactor
detector_width_px = sz_y * oversampleFactor
# We add one extra row and column for the trapezoidal integration method.
x = np.arange(-detector_height_px//2, +detector_height_px//2 + detector_height_px%2 + 1, 1) + 0.
y = np.arange(-detector_width_px//2, +detector_width_px//2 + detector_width_px%2 + 1, 1) + 0.
x *= l_px_m / oversampleFactor
y *= l_px_m / oversampleFactor
Y, X = np.meshgrid(y, x)

I_grid = I(X,Y)
nan_idx = np.where(np.isnan(I_grid))
if nan_idx[0].shape != (0,):
	I_grid[nan_idx[0][0],nan_idx[1][0]] = I_0 # removing the NaN in the centre of the image if necessary

# Numerically integrating using Simpson's (or trapezoidal) rule:
# In x first, then in y.
# First, just integrate over the whole field to see if it returns the expected vale.
res = integrate.cumtrapz(I_grid, dx = l_px_m/oversampleFactor, axis = 0, initial = 0)
res2 = integrate.cumtrapz(res[-1,:], dx = l_px_m/oversampleFactor, initial = 0)
print 'Total count using cumtrapz:', res2[-1]
print 'Total count using I * A_px:',sum(I_grid.flatten())*l_px_m**2/oversampleFactor**2

# Now, we want to repeat this process for each pixel.
photon_grid = np.zeros((sz_x,sz_y))
cumsum = 0
for j in range(sz_y):
	for k in range(sz_x):
		px_grid = I_grid[oversampleFactor*k:oversampleFactor*k+oversampleFactor+1,oversampleFactor*j:oversampleFactor*j+oversampleFactor+1]
		# Checking
		cumsum += sum(px_grid.flatten())*l_px_m**2/oversampleFactor**2
		res1 = integrate.cumtrapz(px_grid, dx = l_px_m/oversampleFactor, axis = 0, initial = 0)
		res2 = integrate.cumtrapz(res1[-1,:], dx = l_px_m/oversampleFactor, initial = 0)
		photon_grid[k,j] = res2[-1]
print 'Total count after oversampling and numerical integration:',sum(photon_grid.flatten())
print 'Total count using sum(px_grid * A_px):', cumsum

# Comparing...
plt.figure(figsize=(3*FIGSIZE,FIGSIZE))
plt.subplot(1,3,1)
plt.imshow(I_grid)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title('Intensity')
plt.subplot(1,3,2)
plt.imshow(photon_grid)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title('Pixel counts')
plt.subplot(1,3,3)
plt.imshow(I_grid*l_px_m**2/oversampleFactor**2)
plt.colorbar(fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
plt.title('Oversampled approximated pixel counts')
plt.show()

count_cumtrapz, I, P_0, P_sum, I_0 = airyDisc(wavelength_m = wavelength_m, f_ratio = f_ratio, l_px_m = l_px_m, detector_size_px = detector.size_px, plotIt = True)