############################################################################################
#
# 	File:		ossim.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	For generating class instances of optical systems that will be used a lot.
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
#
#	This file is part of linguinesim.
#
#	linguinesim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	linguinesim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with linguinesim.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################
from __future__ import division
import linguinesim
from linguinesim.apdsim import *

def anu23mTelescope():
	"""
		Make a TelescopeClass instance corresponding to the ANU 2.3 m.

		TODO: make a telescope using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""
	tel = linguinesim.Telescope(
		efl_m = 6010 / 700 * 4.840, # Taken from Gabe Bloxham's and John Hart's notes 
		T = 273 + 10
		)
	
	# M1
	tel.addMirror(
		R_outer_m = 2.337 / 2,
		R_inner_m = 0.229 / 2
		) # Check that the effective emissivity works correctly.

	# M2 
	# For now, don't worry about the baffle. Talk to Rob first
	tel.addMirror(
		R_outer_m = 0.355 / 2
		)

	# M3
	# We don't really care about the radius of M3 since it never gets used. 
	tel.addMirror(
		R_outer_m = 0.355 / 2
		)

	# For the spider, we pretend that it's a mirror. 
	tel.addSpider(
		A_spider_m2 = 4 * 0.012 * (tel.mirrors[0].R_outer_m - tel.mirrors[1].R_outer_m),	# Total area of spider (m^2)
		eps_spider = AL_EMISSIVITY
		)

	return tel

############################################################################################
def saphiraDetector():
	"""
		Make a DetectorClass instance corresponding to the ANU 2.3 m.

		TODO: make a detector using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""
	return linguinesim.Detector(
		height_px = 256,				# height (pixels)
		width_px = 320,					# width (pixels)
		l_px_m = 24e-6,					# pixel width (m)
		wavelength_cutoff = 2.5e-6,		# cutoff wavelength (m)
		RN = 9,							# ? sqrt(e/pixel) rms
		gain = 50,						# ? avalanche gain
		dark_current = 0.03 * 50,		# ? e/second/pixel; worst-case
		saturation = 2**16 - 1,			# ? detector saturation limit
		adu_gain = 1/2.9,				# electrons per ADU at readout
		qe = 0.9						# quantum efficiency
		)

############################################################################################
def saphiraCryostat():
	print 'TODO: Cryostat temperature needs updating!'
	
	return linguinesim.Cryostat(
		T = 172.372,
		Tr_win = 0.98,
		Omega = np.pi,
		eps_wall = 1.0
		)

############################################################################################
def ssoSky():
	return linguinesim.Sky(
			magnitude_system = 'AB',
			brightness = {
				'J' : 15,
				'H' : 13.7,
				'K' : 12.5
			},
			T = 273,
			eps = getSkyEps()
		)

############################################################################################
def linguineOpticalSystem():
	"""
		Make an OpticalSystemClass instance corresponding to the ANU 2.3 m.

		TODO: make an optical system using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""

	return linguinesim.OpticalSystem(
		telescope = anu23mTelescope(),
		detector = saphiraDetector(),
		cryostat = saphiraCryostat(),
		sky = ssoSky()
		)
