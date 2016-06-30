############################################################################################
#
# 	File:		galsim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	For simulating images of galaxies. 
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
from apdsim import *

def sersic(n, R_e, R,
	# R_trunc = np.inf,
	mu_e = None,
	I_e = None
	):
	" Returns the Sersic profile at radial distance R in units of surface brightness with index n given half-light radius R_e and mu(R=R_e) = mu_e. "

	# Calculating the constraining parameter.
	if mu_e == None and I_e == None:
		print 'ERROR: at least one of mu_e or I_e must be specified!'
		return
	if mu_e == None:
		mu_e = - 2.5 * np.log10(I_e)
	elif I_e == None:
		I_e = np.power(10, - mu_e / 2.5)
	else:
		print 'WARNING: you have specified both mu_e and I_e, so mu and I might not be consistent! Proceed with caution...'

	# Calculating b_n given the Sersic index n.
	if n > 0.5 and n < 8:
		b_n = 1.9992 * n - 0.3271
	elif n >= 8:
		b_n = 2 * n - 1 / 3
	else:
		print 'ERROR: haven\'t implemented n < 0.5 yet! Check again later...'
		return

	I = I_e * np.exp(- b_n * (np.power(R/R_e, 1/n) - 1))
	mu = - 2.5 * np.log10(I)
	return R, mu, I

############################################################################################
def sersic2D(n, R_e,
	theta_rad = 0,		# Angle between major axis and detector horizontal (radians)
	i_rad = 0,			# Inclination angle (radians; face-on corresponds to i = 0)
	R_max = None,		# Plotting limit. Default is 20 * R_e
	R_trunc = np.inf, 	# Disc truncation radius
	gridsize = 500,		# Number of returned grid points
	mu_e = None,
	I_e = None,
	plotIt = False
	):	
	" Returns 2D Sersic intensity and surface brightness plots. "
	if R_max == None:
		R_max = 20 * R_e

	# Making a 2D intensity plot of the galaxy given its inclination and orientation.
	# scaleFactor = np.sqrt(2)*np.cos(np.pi/4 - theta_rad)
	dR = 2 * R_max / 500
	scaleFactor = 2
	imsize = gridsize*scaleFactor
	r = np.linspace(-R_max*scaleFactor, +R_max*scaleFactor, imsize)
	X, Y = np.meshgrid(r, r)
	R = rotateAndCrop(image_in_array = np.sqrt(X * X + Y * Y / (np.cos(i_rad) * np.cos(i_rad))), angle=theta_rad * 180 / np.pi, cropArg=(imsize-gridsize)//2)
	R, mu_map, I_map = sersic(n=n, R_e=R_e, R=R, mu_e=mu_e, I_e=I_e)
	mu_map[R>R_trunc] = np.inf
	I_map[R>R_trunc] = 0

	if plotIt:
		plt.figure(figsize=(2*figsize,figsize))
		plt.subplot(1,2,1)
		plt.imshow(I_map,norm=LogNorm())
		plt.colorbar()
		plt.title('2D intensity map')
		plt.subplot(1,2,2)
		plt.imshow(mu_map)
		plt.colorbar()
		plt.title('2D surface brightness map')
		plt.show()

	return R, dR, I_map, mu_map,