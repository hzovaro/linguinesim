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

	opticalsystem = OpticalSystem(
		telescope = eos18mTelescope(),	
		detector = nuvuDetector(),
		sky = msoSky(),
		plate_scale_as_px = 0.044
		)

	return opticalsystem

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
		height_px = 512,				# height (pixels)
		width_px = 512,					# width (pixels)
		l_px_m = 16e-6,					# pixel width (m)
		wavelength_cutoff = 1.1e-6,		# cutoff wavelength (m)
		RN = 0.1,						# ? sqrt(e/pixel) rms
		gain = 1000,					# ? EM gain
		cic = 0.001,					# clock-induced charge (e/pixel/frame)
		dark_current = 0.0002,			# ? e/second/pixel; worst-case
		saturation = 2**16 - 1,			# ? detector saturation limit
		adu_gain = 1/2.9,				# electrons per ADU at readout
		qe = 0.9,						# quantum efficiency
		fps = 60						# framerate
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
def aoiAoSystem(wave_height_px,
	wavelength_wfs_m = 589e-9,
	wavelength_science_m = 800e-9,
	rng_seed = 1
	):
	try:
		from aosim.pyxao import wavefront, deformable_mirror, wfs, ao_system, atmosphere, seeing_limited_system
	except:
		print("WARNING: I cannot import pyxao - I am returning instead")
		return

	"""
		Make an AO system instance for AOI.
	"""
	wavefrontPupil = {	
		'type':'annulus',
		'dout': 1.752,
		'din' : 0.250
	}

	# Wave parameters
	m_per_px = wavefrontPupil['dout'] / wave_height_px		# Physical mapping of wave onto primary mirror size

	# Imaging parameters

	# AO system parameters
	N_actuators = 17
	N_lenslets = 16
	actuator_pitch_m = wavefrontPupil['dout'] / N_actuators
	lenslet_pitch_m = wavefrontPupil['dout'] / N_lenslets
	ho_loop_rate = 2000

	dm_geometry = 'square'
	wfs_geometry = 'square'
	central_actuator = 'true'
	central_lenslet = 'false'
	
	edge_radius = 1.4	
	influence_fun = 'gaussian'
	pokeStroke = 1e-7	

	# Atmospheric conditions at MSO
	wavelength_ref_m = 550e-9		# Wavelength reference for Fried parameter (Bennet et al. 2012)
	r0_ref_m = 10e-2
	r0_wfs = np.power((wavelength_wfs_m / wavelength_ref_m), 1.2) * r0_ref_m
	r0_science = np.power((wavelength_science_m / wavelength_ref_m), 1.2) * r0_ref_m 

	# v_wind_m = 10			# Turbulent layer wind speed (m/s)
	# wind_angle_deg = 0.0	# Turbulent layer wind direction (rad)
	# elevation_m = 1000		# Turbulent layer elevation (m)
	# airmass = 1.0			# Airmass

	# These values from Bennet et al. 2012
	r0_ref_m = [r0_ref_m * l for l in [4, 2, 3, 1]]
	v_wind_m = [10, 5, 60, 90]				# Turbulent layer wind speed (m/s)
	wind_angle_deg = [0.0, np.pi/6, 0, 0]	# Turbulent layer wind direction (rad)
	elevation_m = [0, 400, 6000, 19000]		# Turbulent layer elevation (m)
	airmass = 1.0			# Airmass

	# Setting up AO system
	wf_wfs = wavefront.Wavefront(wave = wavelength_wfs_m, m_per_px = m_per_px, sz = wave_height_px, pupil = wavefrontPupil)
	wf_science = wavefront.Wavefront(wave = wavelength_science_m, m_per_px = m_per_px, sz = wave_height_px, pupil = wavefrontPupil)
	wavefronts_dm = [wf_wfs, wf_science] 	# Wavefronts corrected by the DM (in a CL AO system, it's all of them!)
	wavefronts_wfs = [wf_wfs]				# Wacefronts sensed by the WFS
	psf_ix = 1		# Index in the list of wavefronts passed to the DM instance corresponding to the PSF to return

	dm = deformable_mirror.DeformableMirror(
		wavefronts = wavefronts_dm, 
		influence_function = 'gaussian', 
		central_actuator = central_actuator, 
		actuator_pitch = actuator_pitch_m, 
		geometry = dm_geometry, 
		edge_radius = 1.4)

	sh_wfs = wfs.ShackHartmann(
		wavefronts = wavefronts_wfs, 
		lenslet_pitch = lenslet_pitch_m, 
		geometry = wfs_geometry, 
		central_lenslet = central_lenslet, 		
		sampling = 1)
		# fratio = wfs_fratio)

	aoi_ao_system = ao_system.SCFeedBackAO(dm = dm, wfs = sh_wfs, image_ixs = psf_ix)
	
	# The atmosphere is a PHASE SCREEN
	atm = atmosphere.Atmosphere(sz = wave_height_px, 
		m_per_px = m_per_px,
		elevations = elevation_m, 
		r_0 = r0_ref_m, 
		wave_ref = wavelength_ref_m, 
		angle_wind = wind_angle_deg,
		v_wind = v_wind_m, 
		airmass = airmass, 
		seed = rng_seed)

	wf_wfs.add_atmosphere(atm)
	wf_science.add_atmosphere(atm)

	# aoi_ao_system.response_matrix = np.load("/Users/azovaro/python/Modules/aosim/pyxao/aoi_response_matrix.npz")['response_matrix']
	# aoi_ao_system.reconstructor = np.load("/Users/azovaro/python/Modules/aosim/pyxao/aoi_reconstructor_matrix.npz")['reconstructor']

	# Either a full path can be given, or else the file is assumed to be located in the same directory from which the script calling this method is being called.
	aoi_ao_system.response_matrix = np.load("aoi_response_matrix.npz")['response_matrix']
	aoi_ao_system.reconstructor = np.load("aoi_reconstructor_matrix.npz")['reconstructor']

	return aoi_ao_system

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
		eps_spider = 1.0
		# eps_spider = AL_EMISSIVITY
		)

	return tel

############################################################################################
def saphiraDetector():
	"""
		Make a DetectorClass instance corresponding to the ANU 2.3 m.

		TODO: make a detector using this function, compare all of its parameters to that made in the sysparams class (to be removed)
	"""
	# print("TODO: fix gain. Avalanche gain currently set to 1.")
	return Detector(
		height_px = 256,				# height (pixels)
		width_px = 320,					# width (pixels)
		l_px_m = 24e-6,					# pixel width (m)
		wavelength_cutoff = 2.5e-6,		# cutoff wavelength (m)
		RN = 9,							# ? sqrt(e/pixel) rms
		gain = 500,					# ? avalanche gain
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
			# Source: GMTIFS On-line Exposure Time Calculator (ETC) and GMTIFSsim input data, rev 1.1
			brightness = {
				'J' : 16.61,
				'H' : 15.49,
				'K' : 14.45	# 15.35/14.45 for winter/summer
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

############################################################################################
def linguineAoSystem(wave_height_px,
	rng_seed = 1
	):
	"""
		Make an AO system instance for the SAPHIRA-2.3 m telescope system. This is only used to generate seeing-limited and diffraction-limited PSFs as the telescope doesn't have an AO system.
	"""
	try:
		from aosim.pyxao import wavefront, deformable_mirror, wfs, ao_system, atmosphere, seeing_limited_system
	except:
		print("WARNING: I cannot import pyxao - I am returning instead")
		return

	wavefront_pupil = {	
		'type':'annulus',
		'dout': 2.337,
		'din' : 0.335,
	}

	# Wave parameters
	m_per_px = wavefront_pupil['dout'] / wave_height_px		# Physical mapping of wave onto primary mirror size

	# Atmospheric conditions at SSO
	wavelength_ref_m = 500e-9		# Wavelength reference for Fried parameter
	r0_ref_m = 10e-2

	v_wind_m = 10			# Turbulent layer wind speed (m/s)
	wind_angle_deg = 0.0	# Turbulent layer wind direction (rad)
	elevation_m = 1000		# Turbulent layer elevation (m)
	airmass = 1.0			# Airmass

	wavelengths_science_m = []
	wavelength_ixs = {}
	r0_science_m = []
	wavefronts = []
	bands = ['J', 'H', 'K']
	for k in range(len(bands)):
		band = bands[k]
		wavelength_ixs[band] = k
		wavelengths_science_m.append(FILTER_BANDS_M[band][0])
		r0_science_m.append(np.power((FILTER_BANDS_M[band][0] / wavelength_ref_m), 1.2) * r0_ref_m)
		# Making a wavefront instance
		wavefronts.append(
			wavefront.Wavefront(wave = FILTER_BANDS_M[band][0],
				m_per_px = m_per_px,
				sz = wave_height_px,
				pupil = wavefront_pupil)
			)	
	
	# The atmosphere is a PHASE SCREEN: not dependent on wavelegnth!
	atm = atmosphere.Atmosphere(sz = wave_height_px, 
		m_per_px = m_per_px,
		elevations = elevation_m, 
		r_0 = r0_ref_m, 
		wave_ref = wavelength_ref_m, 
		angle_wind = wind_angle_deg,
		v_wind = v_wind_m, 
		airmass = airmass, 
		seed = rng_seed)

	for w in wavefronts:
		w.add_atmosphere(atm)

	linguine_ao_system = seeing_limited_system.SeeingLimitedOpticalSystem(wavefronts = wavefronts, wavelength_ixs = wavelength_ixs)

	return linguine_ao_system