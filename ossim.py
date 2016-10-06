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
from __future__ import division, print_function
import numpy as np

from linguineglobals import *
from telescopeclass import Telescope
from detectorclass import Detector
from cryostatclass import Cryostat
from opticalsystemclass import OpticalSystem
from skyclass import Sky
from galaxyclass import Galaxy
import etc

############################################################################################
def aoiOpticalSystem():

	os = OpticalSystem(
		telescope = eos18mTelescope(),
		detector = nuvuDetector(),
		sky = msoSky()
		)
	# Manually setting the plate scale
	os.plate_scale_as_px = 0.044
	os.plate_scale_rad_px = np.deg2rad(os.plate_scale_as_px * 3600)
	os.FoV_height_as = os.detector.height_px * os.plate_scale_as_px
	os.FoV_width_as = os.detector.width_px * os.plate_scale_as_px
	os.FoV_height_rad = os.detector.height_px * os.plate_scale_rad_px
	os.FoV_width_rad = os.detector.width_px * os.plate_scale_rad_px
	os.omega_px_as2 = os.plate_scale_as_px**2
	os.omega_px_sr = os.plate_scale_rad_px**2
	os.etendue = os.omega_px_sr * os.telescope.A_collecting_m2

	return os

############################################################################################
def eos18mTelescope():
	tel = Telescope(
		efl_m = np.inf, 		# For AOI, the plate scale is determined by the optics just before the imager, not by the diameter of the telescope itself (which is afocal anyway)
		T = 273 + 20	# Room temperature (we're in the optical so doesn't matter that much)
		)

	# M1
	tel.addMirror(
		R_outer_m = 1.752 / 2,
		R_inner_m = 0.250 / 2
		)

	# M2 through M7
	for k in range(6):		
		tel.addMirror(R_outer_m = 0.250 / 2)

	return tel

############################################################################################
def nuvuDetector():
	return Detector(
		height_px = 256,				# height (pixels)
		width_px = 256,					# width (pixels)
		l_px_m = 16e-6,					# pixel width (m)
		wavelength_cutoff = 1.1e-6,		# cutoff wavelength (m)
		RN = 0.1,						# ? sqrt(e/pixel) rms
		gain = 1,						# ? avalanche gain
		cic = 0.001,					# clock-induced charge (e/pixel/frame)
		dark_current = 0.0002,			# ? e/second/pixel; worst-case
		saturation = 2**16 - 1,			# ? detector saturation limit
		adu_gain = 1/2.9,				# electrons per ADU at readout
		qe = 0.9						# quantum efficiency
		)

############################################################################################
def msoSky():
	return Sky(
		magnitude_system = 'AB',
			brightness = {
				'J' : 15,
				'H' : 13.7,
				'K' : 12.5
			},
			T = 273,
			eps = etc.getSkyEps()
		)

############################################################################################
def aoiAoSystem():
	"""
		Make an AO system instance for AOI.
	"""
	

############################################################################################
############################################################################################
def anu23mTelescope():
	"""
		Make a TelescopeClass instance corresponding to the ANU 2.3 m.

		TODO: make a telescope using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""
	tel = Telescope(
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
		R_outer_m = 0.355 / 2,
		eps_eff = 1.0
		)

	# M3
	# We don't really care about the radius of M3 since it never gets used. 
	tel.addMirror(
		R_outer_m = 0.355 / 2
		)

	# For the spider, we pretend that it's a mirror. 
	tel.addSpider(
		A_spider_m2 = 4 * 0.012 * (tel.mirrors[0].R_outer_m - tel.mirrors[1].R_outer_m),	# Total area of spider (m^2)
		# eps_spider = 1.0
		eps_spider = AL_EMISSIVITY
		)

	return tel

############################################################################################
def saphiraDetector():
	"""
		Make a DetectorClass instance corresponding to the ANU 2.3 m.

		TODO: make a detector using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""
	print("TODO: fix gain. Avalanche gain currently set to 1.")
	return Detector(
		height_px = 256,				# height (pixels)
		width_px = 320,					# width (pixels)
		l_px_m = 24e-6,					# pixel width (m)
		wavelength_cutoff = 2.5e-6,		# cutoff wavelength (m)
		RN = 9,							# ? sqrt(e/pixel) rms
		gain = 1,						# ? avalanche gain
		dark_current = 0.03,			# ? MULTIPLY BY GAIN!! e/second/pixel; worst-case
		saturation = 2**16 - 1,			# ? detector saturation limit
		adu_gain = 1/2.9,				# electrons per ADU at readout
		qe = 0.9						# quantum efficiency
		)

############################################################################################
def saphiraCryostat():
	print("TODO: Cryostat temperature needs updating!")
	
	return Cryostat(
		T = 172.372,
		Tr_win = 0.98,
		Omega = np.pi,
		eps_wall = 1.0
		)

############################################################################################
def ssoSky():
	return Sky(
			magnitude_system = 'AB',
			brightness = {
				'J' : 15,
				'H' : 13.7,
				'K' : 12.5
			},
			T = 273,
			eps = etc.getSkyEps()
		)

############################################################################################
def linguineOpticalSystem():
	"""
		Make an OpticalSystemClass instance corresponding to the ANU 2.3 m.

		TODO: make an optical system using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""

	return OpticalSystem(
		telescope = anu23mTelescope(),
		detector = saphiraDetector(),
		cryostat = saphiraCryostat(),
		sky = ssoSky()
		)
