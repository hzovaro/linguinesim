####################################################################################################
#
# 	File:		etc.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	An exposure time calculator (ETC) for a telescope-detector system in the near-infrared.
#	The detector is assumed to be housed in a cryostat. 
#	
#	The returned SNR is per pixel (for now).
#	
#	The user specifies the exposure time, the surface brightness (in Vega or AB magnitudes)
#	and the imaging band (J, H or K).
#
#	Copyright (C) 2016 Anna Zovaro
#
####################################################################################################
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
####################################################################################################
from __future__ import division, print_function
import miscutils as mu
import numpy as np
import pdb
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib import rc
rc('image', interpolation='none', cmap = 'binary_r')

import scipy.constants

import json	

from linguineglobals import *
import etcutils
####################################################################################################
def exposureTimeCalc(band, t_exp, optical_system,
		surface_brightness = None,
		magnitude_system = None,
		printIt = True
	):
	"""
		An exposure time calculator. 
		Outputs are expressed per pixel.
		Noise terms include contributions from
		- the telescope (spider structure and mirrors, radiating as blackbodies at telescope.T)
		- the sky (radiating as a blackbody at sky.T)
		- the cryostat (radiating as a blackbody at cryostat.T)
		- detector read noise 
		- detector dark current.

		In the J and H bands, empirical sky brightness values from the sky module are used instead of theoretical estimates in the SNR result.
		In the K band, the sky and telescope values are used instead.
		In either case, all noise counts (empirical and theoretical) are returned in their respective fields in the output dictionary.

		The returned counts/sec 

	"""
	detector = optical_system.detector
	telescope = optical_system.telescope
	cryostat = optical_system.cryostat
	sky = optical_system.sky

	wavelength_eff = FILTER_BANDS_M[band][0]
	bandwidth = FILTER_BANDS_M[band][1]
	wavelength_min = FILTER_BANDS_M[band][2]
	wavelength_max = FILTER_BANDS_M[band][3]

	################################################################################################
	# Noise terms in the SNR calculation
	################################################################################################
	
	""" Signal photon flux """
	# Given the surface brightness of the object, calculate the source electrons/s/pixel.
	if surface_brightness != None:
		if magnitude_system != None:
			# Here, if the input is given in mag/arcsec^2, then we need Sigma_source_e to be returned in units of electrons/s/pixel. 
			Sigma_source_e = etcutils.surface_brightness2countRate(mu = surface_brightness, 
				wavelength_m = wavelength_eff, 
				bandwidth_m = bandwidth, 
				plate_scale_as_px = optical_system.plate_scale_as_px, 
				A_tel = telescope.A_collecting_m2, 
				tau = telescope.tau * cryostat.Tr_win,
				qe = detector.qe,
				gain = 1,
				magnitude_system = magnitude_system
			)
		else:
			print('ERROR: you must specify a magnitude system for the source!')
			return
	else:
		Sigma_source_e = 0

	""" Cryostat photon flux """
	Sigma_cryo = getCryostatTE(optical_system=optical_system)

	""" Telescope thermal background photon flux """
	Sigma_tel = getTelescopeTE(optical_system=optical_system, plotit=False)[band]

	""" Sky thermal background photon flux """
	Sigma_sky_thermal = getSkyTE(optical_system=optical_system, plotit=False)[band]

	""" Empirical sky background flux """
	Sigma_sky_emp = etcutils.surface_brightness2countRate(mu = sky.brightness[band], 
		wavelength_m = wavelength_eff, 
		bandwidth_m = bandwidth, 
		plate_scale_as_px = optical_system.plate_scale_as_px, 
		A_tel = telescope.A_collecting_m2, 
		tau = telescope.tau * cryostat.Tr_win,
		qe = detector.qe,
		gain = 1,
		magnitude_system = sky.magnitude_system
	)

	""" Total sky background """
	if band == 'K':
		# In the K band, thermal emission from the sky is dominated by the telescope and sky thermal emission
		Sigma_sky = Sigma_sky_thermal + Sigma_tel
	else:
		# In the J and H bands, OH emission dominates; hence empirical sky brightness values are used instead.
		Sigma_sky = Sigma_sky_emp

	""" Dark current """
	# Be careful about gain compensation! 
	Sigma_dark = detector.dark_current

	################################################################################################
	# Calculating the SNR
	################################################################################################
	N_source = Sigma_source_e * t_exp
	N_dark = Sigma_dark * t_exp
	N_cryo = Sigma_cryo * t_exp
	N_sky = Sigma_sky * t_exp
	N_tel = Sigma_tel * t_exp
	N_sky_emp = Sigma_sky_emp * t_exp
	N_sky_thermal = Sigma_sky_thermal * t_exp
	N_RN = detector.RN**2

	SNR_unity_gain = N_source / np.sqrt(N_source + N_dark + N_cryo + N_sky + N_RN)
	SNR_gain_multiplied = N_source * detector.gain / np.sqrt(N_source* detector.gain + N_dark* detector.gain + N_cryo* detector.gain + N_sky* detector.gain + N_RN)
	################################################################################################

	etc_output = {
		# Input parameters
		't_exp' : t_exp,
		'band' : band,
		'surface_brightness' : surface_brightness,
		'magnitude_system' : magnitude_system,		
		
		'gain-multiplied' : {
			# Noise standard deviations PER PIXEL
			# Poisson distribution, so the nosie scales as the square root of the total number of photons
			'sigma_source' : np.sqrt(N_source * detector.gain),
			'sigma_dark' : np.sqrt(N_dark * detector.gain),
			'sigma_cryo' : np.sqrt(N_cryo * detector.gain),
			'sigma_sky' : np.sqrt(N_sky * detector.gain),
			'sigma_sky_emp' : np.sqrt(N_sky_emp * detector.gain),
			'sigma_tel' : np.sqrt(N_tel * detector.gain),
			'sigma_sky_thermal' : np.sqrt(N_sky_thermal * detector.gain),
			'sigma_RN' : detector.RN,

			# Total electron count PER PIXEL at the given exposure time WITH gain multiplication
			'N_source' : N_source * detector.gain,
			'N_dark' : N_dark * detector.gain,
			'N_cryo' : N_cryo * detector.gain,
			'N_sky_emp' : N_sky_emp * detector.gain,
			'N_sky_thermal' : N_sky_thermal * detector.gain,
			'N_sky' : N_sky * detector.gain,
			'N_tel' : N_tel * detector.gain,
			'N_RN' : N_RN,
			'SNR' : SNR_gain_multiplied
		},

		# Total electron count PER PIXEL with NO GAIN MULTIPLICATION
		'unity gain' : {
			# Noise standard deviations PER PIXEL
			# Poisson distribution, so the nosie scales as the square root of the total number of photons
			'sigma_source' : np.sqrt(N_source),
			'sigma_dark' : np.sqrt(N_dark),
			'sigma_cryo' : np.sqrt(N_cryo),
			'sigma_sky' : np.sqrt(N_sky),
			'sigma_sky_emp' : np.sqrt(N_sky_emp),
			'sigma_tel' : np.sqrt(N_tel),
			'sigma_sky_thermal' : np.sqrt(N_sky_thermal),
			'sigma_RN' : detector.RN,

			# Total electron count PER PIXEL at the given exposure time WITHOUT gain multiplication
			'N_source' : N_source,
			'N_dark' : N_dark,
			'N_cryo' : N_cryo,
			'N_sky_emp' : N_sky_emp,
			'N_sky_thermal' : N_sky_thermal,
			'N_sky' : N_sky,
			'N_tel' : N_tel,
			'N_RN' : N_RN,
			'SNR' : SNR_unity_gain
		}
	}

	if printIt:
		print(json.dumps(etc_output, indent=4, sort_keys=True))
		mu.println()
		if band == "K":
			print("WARNING: Imaging in {:}-band: sky background includes modelled thermal emission contributions from sky and telescope".format(band))
		else:
			print("WARNING: Imaging in {:}-band: sky background taken from empirical sky brightness measurements".format(band))
		print("Origin\t\tExpected count per pixel in frame (unity gain), t_exp = {:.2f} s".format(t_exp))
		mu.println()
		print("Source\t\t{:10.10g}".format(N_source))
		print("Dark current\t{:10.10g}".format(N_dark))
		print("Cryostat\t{:10.10g}".format(N_cryo))
		print("Sky\t\t{:10.10g}".format(N_sky))
		print("Read noise\t{:10.10g}".format(N_RN))
		mu.println()
		print("SNR (gain-multiplied)\t\t{:10.10g}".format(SNR_gain_multiplied))
		print("SNR (unity gain)\t\t{:10.10g}".format(SNR_unity_gain))
		mu.println()

	return etc_output

####################################################################################################
def getCryostatTE(optical_system):
	"""
		Return the expected counts/second/pixel on a detector resulting from thermal emission from the cryostat walls.

		Inputs:
		------------
		optical_system: OpticalSystemClass instance 
			An OpticalSystemClass with a CryostatClass member.

		Returns:
		------------
		I_cryo: float
			The expected count (in electrons/second/pixel) in the detector resulting from the thermal emission from the cryostat walls.
	"""
	detector = optical_system.detector
	cryostat = optical_system.cryostat

	return etcutils.thermalEmissionIntensity(
			T = cryostat.T,
			A = detector.A_px_m2,
			wavelength_min = 0.0,
			wavelength_max = detector.wavelength_cutoff,
			Omega = cryostat.Omega,
			eps = cryostat.eps_wall,
			eta = detector.qe
			)



####################################################################################################

def getSkyTE(optical_system,
	plotit=True):

	" Sky thermal background photon flux in the J, H and K bands "

	detector = optical_system.detector
	telescope = optical_system.telescope
	cryostat = optical_system.cryostat
	sky = optical_system.sky

	I_sky = {
		'J' : 0.0,
		'H' : 0.0,
		'K' : 0.0
	}

	# Atmospheric properties	
	# eps_sky = getSkyEps()

	for key in I_sky:
		wavelength_min = FILTER_BANDS_M[key][2]
		wavelength_max = FILTER_BANDS_M[key][3]
		I_sky[key] = etcutils.thermalEmissionIntensity(
			T = sky.T, 
			wavelength_min = wavelength_min, 
			wavelength_max = wavelength_max, 
			Omega = optical_system.omega_px_sr, 
			A = telescope.A_collecting_m2, 
			eps = sky.eps,
			eta = detector.qe * telescope.tau * cryostat.Tr_win
			)

	if plotit:
		D = np.ones(1000)*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		mu.newfigure(1,1)
		plt.plot(wavelengths*1e6, D, 'g--', label=r'Dark current')
		# Imager mode
		for key in I_sky:
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, I_sky[key], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o')
			plt.text(FILTER_BANDS_M[key][0]*1e6, I_sky[key]*5, key)

		plt.yscale('log')
		plt.axis('tight')
		plt.ylim(ymax=100*I_sky['K'],ymin=1e-5)
		plt.legend(loc='lower right')
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel(r'Count ($e^{-}$ s$^{-1}$ pixel$^{-1}$)')
		plt.title(r'Estimated count from sky thermal emission')
		plt.show()

	return I_sky

####################################################################################################
def getTelescopeTE(optical_system,
	plotit=True):

	detector = optical_system.detector
	telescope = optical_system.telescope
	cryostat = optical_system.cryostat
	sky = optical_system.sky

	I_tel = {
		'J' : 0.0,
		'H' : 0.0,
		'K' : 0.0
	}

	for key in I_tel:
		wavelength_min = FILTER_BANDS_M[key][2]
		wavelength_max = FILTER_BANDS_M[key][3]
			
		# Mirrors
		# Assumptions:
		#	1. The area we use for the etendue is the collecting (i.e. reflective) area of the telescope, not the total area.
		#	2. For now we are ignoring the baffle on M2.
		#	3. We are not assuming the worst case for the spider (i.e. it is still substantially reflective). But you should see how substantial of a difference it makes. Always lean towards the worst-case. 
		I_mirrors = 0
		for mirror in telescope.mirrors:
			I_mirrors += etcutils.thermalEmissionIntensity(
				T = telescope.T, 
				wavelength_min = wavelength_min, 
				wavelength_max = wavelength_max, 
				Omega = optical_system.omega_px_sr, 
				A = telescope.A_collecting_m2, 
				eps = mirror.eps_eff,
				eta = detector.qe * cryostat.Tr_win
				)
		
		# Spider 
		if telescope.hasSpider:
			I_spider = etcutils.thermalEmissionIntensity(
					T = telescope.T, 
					wavelength_min = wavelength_min, 
					wavelength_max = wavelength_max, 
					Omega = optical_system.omega_px_sr, 
					A = telescope.A_collecting_m2, 
					eps = telescope.eps_spider_eff,
					eta = telescope.tau * detector.qe * cryostat.Tr_win)\
			  + etcutils.thermalEmissionIntensity(
			  		T = sky.T, 	
			  		wavelength_min = wavelength_min, 
			  		wavelength_max = wavelength_max, 
			  		Omega = optical_system.omega_px_sr, 
			  		A = telescope.A_collecting_m2, 
			  		eps = lambda wavelength_m : (1 - telescope.eps_spider_eff) * sky.eps(wavelength_m),
			  		eta = telescope.tau * detector.qe * cryostat.Tr_win)
		
		# Cryostat window 
		I_window = etcutils.thermalEmissionIntensity(
			T = cryostat.T, 
			wavelength_min = wavelength_min, 
			wavelength_max = wavelength_max, 
			Omega = optical_system.omega_px_sr, 
			A = telescope.A_collecting_m2, 
			eps = cryostat.eps_win,
			eta = detector.qe	# No cryostat window or telescope throughput terms because the radiation from the walls doesn't pass through it
			)

		I_tel[key] = I_mirrors + I_spider + I_window

	if plotit:
		D = np.ones(1000)*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		mu.newfigure(1,1)
		plt.plot(wavelengths*1e6, D, 'g--', label=r'Dark current')

		for key in I_tel:
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, I_tel[key], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o')
			plt.text(FILTER_BANDS_M[key][0]*1e6, I_tel[key]*5, key)

		plt.yscale('log')
		plt.axis('tight')
		plt.ylim(ymax=100*I_tel['K'],ymin=1e-5)
		plt.legend(loc='lower right')
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel(r'Count ($e^{-}$ s$^{-1}$ pixel$^{-1}$)')
		plt.title(r'Estimated count from telescope thermal emission')
		plt.show()

	return I_tel

########################################################################################
def plotBackgroundNoiseSources(optical_system):
	" Plot the empirical sky brightness, thermal sky emission, thermal telescope emission and dark current as a function of wavelength_m "

	detector = optical_system.detector
	telescope = optical_system.telescope
	cryostat = optical_system.cryostat
	sky = optical_system.sky

	counts = {
		'H' : 0,
		'J' : 0,
		'K' : 0
	}
	counts['H'] = exposureTimeCalc(band='H', t_exp=1, optical_system=optical_system)
	counts['J'] = exposureTimeCalc(band='J', t_exp=1, optical_system=optical_system)
	counts['K'] = exposureTimeCalc(band='K', t_exp=1, optical_system=optical_system)
	D = np.ones(1000)*detector.dark_current
	wavelengths = np.linspace(1.0, 2.5, 1000)*1e-6

	# Plotting
	mu.newfigure(1.5,1.5)
	plt.plot(wavelengths*1e6, D, 'g--', label=r'Dark current')
	plotColors = {
		'H' : 'orangered',
		'J' : 'darkorange',
		'K' : 'darkred'
	}
	for key in counts:
		if key == 'J':
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_sky_emp'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o', ecolor=plotColors[key], mfc=plotColors[key], label='Empirical sky background')
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='^', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal sky background')
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_tel'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal telescope background')
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_tel'] + counts[key]['gain-multiplied']['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='x', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal telescope + sky background')
		else:
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_sky_emp'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o', ecolor=plotColors[key], mfc=plotColors[key])
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='^', ecolor=plotColors[key], mfc=plotColors[key])
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_tel'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key])
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_tel'] + counts[key]['gain-multiplied']['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='x', ecolor=plotColors[key], mfc=plotColors[key])

		plt.text(FILTER_BANDS_M[key][0]*1e6, counts[key]['gain-multiplied']['N_sky_emp']*5, key)

	plt.yscale('log')
	plt.axis('tight')
	plt.ylim(ymax=100*counts['K']['gain-multiplied']['N_tel'],ymin=1e-5)
	plt.legend(loc='lower right')
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'Count ($e^{-}$ s$^{-1}$ pixel$^{-1}$)')
	plt.title(r'Expected background noise levels (gain-multiplied by %d)' % detector.gain)
	plt.show()

####################################################################################################
def findCryostatTemp(optical_system, plotit=True):
	"""
		Determine what temperature the cryostat given in the optical system must be so that the detector counts resulting from the thermal emission from the cryostat walls is equivalent to the dark current.

		Inputs:
		------------
		optical_system: OpticalSystemClass instance 
			An OpticalSystemClass with a CryostatClass member.

		Returns:
		------------
		I_cryo: float
			The expected count (in electrons/second/pixel) in the detector resulting from the thermal emission from the cryostat walls.

	"""
	detector = optical_system.detector
	cryostat = optical_system.cryostat

	T_cryo = np.linspace(80, 200, 1000)
	I_cryo = np.zeros(len(T_cryo))
	I_cryo_h = np.zeros(len(T_cryo))

	# For finding the crossover point
	minVal = np.inf	
	minVal_min = np.inf

	# IMPORTANT NOTE: we do NOT multiply by the gain here because (1) it is not yet certain what gain values we will use and (2) the cryostat emission and the dark current are (to first order) both affected by the gain in the same way, so it doesn't matter whether or not we apply the gain here or not AS LONG AS the dark current value stored in the detector instance is the PRE-GAIN value!
	for k in range(T_cryo.size):
		I_cryo[k] = etcutils.thermalEmissionIntensity(
			T = T_cryo[k],
			A = detector.A_px_m2,
			wavelength_min = 0.0,
			wavelength_max = detector.wavelength_cutoff,
			Omega = cryostat.Omega,
			eps = cryostat.eps_wall,
			eta = detector.qe
			)
		# Find the crossing point
		if abs(I_cryo[k] - detector.dark_current) < minVal:
			minVal = abs(I_cryo[k] - detector.dark_current)
			idx = k

		# Worst case: increased cutoff wavelength.
		I_cryo_h[k] = etcutils.thermalEmissionIntensity(
			T = T_cryo[k],
			A = detector.A_px_m2,
			wavelength_min = 0.0,
			wavelength_max = detector.wavelength_cutoff_h,
			Omega = cryostat.Omega,
			eps = cryostat.eps_wall,
			eta = detector.qe
			)
		# Find the crossing point
		if abs(I_cryo_h[k] - detector.dark_current*0.1) < minVal_min:
			minVal_min = abs(I_cryo_h[k] - detector.dark_current*0.1)
			idx_min = k

		print("REMINDER: the detector dark current is currently set to %.4f. Make sure that the stored dark current value is BEFORE gain multiplication or these results are invalid.")

	if plotit:
		D = np.ones(len(T_cryo))*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		mu.newfigure(1.5,1.5)
		plt.rc('text', usetex=True)
		plt.plot(T_cryo, I_cryo, 'r', label='Cryostat thermal emission, $\lambda_c = %.1f \mu$m' % (detector.wavelength_cutoff*1e6))
		plt.plot(T_cryo, I_cryo_h, 'r--', label='Cryostat thermal emission, $\lambda_c = %.1f \mu$m' % (detector.wavelength_cutoff_h*1e6))
		plt.plot(T_cryo, D, 'g', label=r'Dark current')
		plt.plot(T_cryo, D*0.1, 'g--', label=r'Dark current (10\%)')
		plt.plot([T_cryo[idx],T_cryo[idx]], [0, np.max(I_cryo_h)], 'k', label='$T_c = %.3f$ K' % (T_cryo[idx]))
		plt.plot([T_cryo[idx_min],T_cryo[idx_min]], [0, np.max(I_cryo_h)], 'k--', label='$T_c = %.3f$ K (worst-case)' % T_cryo[idx_min])
		plt.yscale('log')
		plt.axis('tight')
		plt.legend(loc='lower right')
		plt.xlabel(r'Temperature (K)')
		plt.ylabel(r'Count ($e^{-}$ s$^{-1}$ pixel$^{-1}$)')
		plt.title('Estimated count from cryostat thermal emission')
	
	return T[idx]

########################################################################################
def getSkyEps():
	# Atmospheric properties	
	fname = 'cptrans_zm_23_10.dat'
	this_dir, this_filename = os.path.split(__file__)
	DATA_PATH = os.path.join(this_dir, 'skytransdata', fname)
	f = open(DATA_PATH, 'r')

	wavelengths_sky = [];
	Tr_sky = [];
	for line in f:
		cols = line.split()
		wavelengths_sky.append(float(cols[0]))
		Tr_sky.append(float(cols[1]))
	Tr_sky = np.asarray(Tr_sky)
	wavelengths_sky = np.asarray(wavelengths_sky) * 1e-6
	eps_sky = lambda wavelength_m: np.interp(wavelength_m, wavelengths_sky, 1 - Tr_sky)
	f.close()

	return eps_sky
