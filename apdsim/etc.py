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
####################################################################################################
from __future__ import division
from apdsim import *

def exposureTimeCalc(band, t_exp, worstCaseSpider,
		surfaceBrightness = None,
		magnitudeSystem = None,
	):
	"""
		An exposure time calculator. 
		Outputs are expressed per pixel.
		Noise terms include contributions from
		- the telescope (spider structure and mirrors, radiating as blackbodies at telescope.T)
		- the sky (radiating as a blackbody at sky.T)
		- the cryostat (radiating as a blackbody at cryo.T)
		- detector read noise 
		- detector dark current.

		In the J and H bands, empirical sky brightness values from the sky module are used instead of theoretical estimates in the SNR result.
		In the K band, the sky and telescope values are used instead.
		In either case, all noise counts (empirical and theoretical) are returned in their respective fields in the output dictionary.

	"""
	wavelength_eff = FILTER_BANDS_M[band][0]
	bandwidth = FILTER_BANDS_M[band][1]
	wavelength_min = FILTER_BANDS_M[band][2]
	wavelength_max = FILTER_BANDS_M[band][3]

	################################################################################################
	# Noise terms in the SNR calculation
	################################################################################################
	
	""" Signal photon flux """
	# Given the surface brightness of the object, calculate the source electrons/s/pixel.
	if surfaceBrightness != None:
		if magnitudeSystem != None:
			# Here, if the input is given in mag/arcsec^2, then we need Sigma_source_e to be returned in units of electrons/s/pixel. 
			Sigma_source_e = surfaceBrightness2countRate(mu = surfaceBrightness, 
				wavelength_m = wavelength_eff, 
				bandwidth_m = bandwidth, 
				plate_scale_as_px = SYSTEM_PLATE_SCALE_AS_PX, 
				A_tel = telescope.A_collecting, 
				tau = telescope.tau * cryo.Tr_win,
				qe = detector.qe,
				gain = detector.gain,
				magnitudeSystem = magnitudeSystem
			)
		else:
			print 'ERROR: you must specify a magnitude system for the source!'
			return
	else:
		Sigma_source_e = 0

	""" Cryostat photon flux """
	Sigma_cryo = getCryostatTE(plotIt=False)

	""" Telescope thermal background photon flux """
	Sigma_tel = getTelescopeTE(T_sky=sky.T, plotIt=False, worstCaseSpider=worstCaseSpider)[band]

	""" Sky thermal background photon flux """
	Sigma_sky_thermal = getSkyTE(T_sky=sky.T, plotIt=False)[band]

	""" Empirical sky background flux """
	Sigma_sky_emp = surfaceBrightness2countRate(mu = sky.brightness[band], 
		wavelength_m = wavelength_eff, 
		bandwidth_m = bandwidth, 
		plate_scale_as_px = SYSTEM_PLATE_SCALE_AS_PX, 
		A_tel = telescope.A_collecting, 
		tau = telescope.tau * cryo.Tr_win,
		qe = detector.qe,
		gain = detector.gain,
		magnitudeSystem = sky.magnitude_system
	)

	""" Total sky background """
	if band == 'K':
		# In the K band, thermal emission from the sky is dominated by the telescope and sky thermal emission
		Sigma_sky = Sigma_sky_thermal + Sigma_tel
	else:
		# In the J and H bands, OH emission dominates; hence empirical sky brightness values are used instead.
		Sigma_sky = Sigma_sky_emp

	""" Dark current """
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
	N_RN = detector.read_noise * detector.read_noise

	SNR = N_source / np.sqrt(N_source + N_dark + N_cryo + N_sky + N_RN)
	################################################################################################

	etc_output = {
		# Input parameters
		't_exp' : t_exp,
		'band' : band,
		'surfaceBrightness' : surfaceBrightness,
		'magnitudeSystem' : magnitudeSystem,
		# Noise standard deviations PER PIXEL
		# Poisson distribution, so the nosie scales as the square root of the total number of photons
		'sigma_source' : np.sqrt(N_source),
		'sigma_dark' : np.sqrt(N_dark),
		'sigma_cryo' : np.sqrt(N_cryo),
		'sigma_sky' : np.sqrt(N_sky),
		'sigma_sky_emp' : np.sqrt(N_sky_emp),
		'sigma_tel' : np.sqrt(N_tel),
		'sigma_sky_thermal' : np.sqrt(N_sky_thermal),
		'sigma_RN' : detector.read_noise,
		# Total electron count PER PIXEL
		'N_source' : N_source,
		'N_dark' : N_dark,
		'N_cryo' : N_cryo,
		'N_sky_emp' : N_sky_emp,
		'N_sky_thermal' : N_sky_thermal,
		'N_sky' : N_sky,
		'N_tel' : N_tel,
		'N_RN' : N_RN,
		# SNR
		'SNR' : SNR
	}

	return etc_output

####################################################################################################
def getCryostatTE(plotIt=True):
	T_cryo = np.linspace(80, 200, 1000)
	I_cryo = np.zeros(len(T_cryo))
	I_cryo_h = np.zeros(len(T_cryo))

	# For finding the crossover point
	minVal = np.inf	
	minVal_min = np.inf

	for k in range(T_cryo.size):
		I_cryo[k] = thermalEmissionIntensity(
			T = T_cryo[k],
			A = detector.A_px_m2,
			wavelength_min = 0.0,
			wavelength_max = detector.wavelength_cutoff,
			Omega = cryo.Omega,
			eps = cryo.eps_wall,
			eta = detector.gain * detector.qe
			)
		I_cryo_h[k] = thermalEmissionIntensity(
			T = T_cryo[k],
			A = detector.A_px_m2,
			wavelength_min = 0.0,
			wavelength_max = detector.wavelength_cutoff_h,
			Omega = cryo.Omega,
			eps = cryo.eps_wall,
			eta = detector.gain * detector.qe
			)
		# Find the crossing point
		if abs(I_cryo[k] - detector.dark_current) < minVal:
			minVal = abs(I_cryo[k] - detector.dark_current)
			idx = k
		# Absolute worst-case
		if abs(I_cryo_h[k] - detector.dark_current*0.1) < minVal_min:
			minVal_min = abs(I_cryo_h[k] - detector.dark_current*0.1)
			idx_min = k

	if plotIt:
		D = np.ones(len(T_cryo))*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		plt.figure(figsize=(FIGSIZE*1.25,FIGSIZE))
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
	
	return I_cryo[idx]

####################################################################################################

def getSkyTE(T_sky, plotIt=True):
	" Sky thermal background photon flux in the J, H and K bands "
	I_sky = {
		'J' : 0.0,
		'H' : 0.0,
		'K' : 0.0
	}

	# Atmospheric properties	
	eps_sky = getSkyEps()

	for key in I_sky:
		wavelength_min = FILTER_BANDS_M[key][2]
		wavelength_max = FILTER_BANDS_M[key][3]
		I_sky[key] = thermalEmissionIntensity(
			T = T_sky, 
			wavelength_min = wavelength_min, 
			wavelength_max = wavelength_max, 
			Omega = OMEGA_PX_RAD, 
			A = telescope.A_collecting, 
			eps = eps_sky
			)
		# Multiply by the gain, QE and telescope transmission to get units of electrons/s/px.
		I_sky[key] = detector.gain * detector.qe * telescope.tau * cryo.Tr_win * I_sky[key]

	if plotIt:
		D = np.ones(1000)*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		plt.figure(figsize=(FIGSIZE,FIGSIZE))
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

def getTelescopeTE(T_sky, plotIt=True, worstCaseSpider=False):
	I_tel = {
		'J' : 0.0,
		'H' : 0.0,
		'K' : 0.0
	}

	for key in I_tel:
		wavelength_min = FILTER_BANDS_M[key][2]
		wavelength_max = FILTER_BANDS_M[key][3]
		
		I_M1 = thermalEmissionIntensity(T = telescope.T, wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = telescope.eps_M1_eff)
		I_M2 = thermalEmissionIntensity(T = telescope.T, wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M2_total_eff, eps = telescope.eps_M2_eff)
		I_M3 = thermalEmissionIntensity(T = telescope.T, wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = telescope.eps_M1_eff)
		
		# Spider (acting as a grey body at both the sky temperature and telescope temperature)
		if worstCaseSpider == False:
			eps_sky = getSkyEps()
			# BEST CASE: assume the spider has a fresh aluminium coating - so 9.1% emissive at sky temp and 90.1% emissive at telescope temp
			I_spider = \
				thermalEmissionIntensity(T = telescope.T, 	wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = telescope.eps_spider_eff)\
			  + thermalEmissionIntensity(T = T_sky, 		wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = lambda wavelength_m : (1 - telescope.eps_spider_eff) * eps_sky(wavelength_m))
		else:
			# WORST CASE: assume the spider is 100% emissive at telescope temperature (i.e. not reflective at all)
			I_spider = \
				thermalEmissionIntensity(T = telescope.T, 	wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = 1.0)\
		
		# Cryostat window (acting as a grey body)
		I_window = thermalEmissionIntensity(T = cryo.T, wavelength_min = wavelength_min, wavelength_max = wavelength_max, Omega = OMEGA_PX_RAD, A = telescope.A_M1_total, eps = cryo.eps_win)

		# Multiply by the gain and QE to get units of electrons/s/px.
		# We don't multiply by the transmission because the mirrors themselves are emitting.
		I_tel[key] = detector.gain * detector.qe * cryo.Tr_win * (I_M1 + I_M2 + I_M3 + I_spider) + detector.gain * detector.qe * I_window

	if plotIt:
		D = np.ones(1000)*detector.dark_current
		wavelengths = np.linspace(0.80, 2.5, 1000)*1e-6

		# Plotting
		plt.figure(figsize=(FIGSIZE,FIGSIZE))
		plt.plot(wavelengths*1e6, D, 'g--', label=r'Dark current')

		for key in FILTER_BANDS_M:
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
def plotBackgroundNoiseSources():
	" Plot the empirical sky brightness, thermal sky emission, thermal telescope emission and dark current as a function of wavelength_m "
	counts = {
		'H' : 0,
		'J' : 0,
		'K' : 0
	}
	counts['H'] = exposureTimeCalc(band='H', t_exp=1)
	counts['J'] = exposureTimeCalc(band='J', t_exp=1)
	counts['K'] = exposureTimeCalc(band='K', t_exp=1)
	N_tel_worstcase = getTelescopeTE(worstCaseSpider=True, plotIt=False)
	D = np.ones(1000)*detector.dark_current
	wavelengths = np.linspace(1.0, 2.5, 1000)*1e-6

	# Plotting
	plt.figure(figsize=(FIGSIZE,FIGSIZE))
	plt.plot(wavelengths*1e6, D, 'g--', label=r'Dark current')
	plotColors = {
		'H' : 'orangered',
		'J' : 'darkorange',
		'K' : 'darkred'
	}
	for key in counts:
		if key == 'J':
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_sky_emp'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o', ecolor=plotColors[key], mfc=plotColors[key], label='Empirical sky background')
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='^', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal sky background')
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_tel'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal telescope background')
			# eb=plt.errorbar(FILTER_BANDS_M[key][0]*1e6, N_tel_worstcase[key], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key], label='Thermal telescope background (worst-case)')
			# eb[-1][0].set_linestyle('--')
		else:
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_sky_emp'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='o', ecolor=plotColors[key], mfc=plotColors[key])
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_sky_thermal'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='^', ecolor=plotColors[key], mfc=plotColors[key])
			plt.errorbar(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_tel'], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key])
			# eb=plt.errorbar(FILTER_BANDS_M[key][0]*1e6, N_tel_worstcase[key], 0, FILTER_BANDS_M[key][1]/2*1e6, fmt='*', ecolor=plotColors[key], mfc=plotColors[key])
			# eb[-1][0].set_linestyle('--')

		plt.text(FILTER_BANDS_M[key][0]*1e6, counts[key]['N_sky_emp']*5, key)

	plt.yscale('log')
	plt.axis('tight')
	plt.ylim(ymax=100*counts['K']['N_tel'],ymin=1e-5)
	plt.legend(loc='lower right')
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'Count ($e^{-}$ s$^{-1}$ pixel$^{-1}$)')
	plt.title(r'Expected background noise levels')
	plt.show()

########################################################################################
def thermalEmissionIntensity(		
	T,					# Emission source temperature
	wavelength_min,		# Integrand interval lower limit
	wavelength_max,		# Integrand interval upper limit
	Omega,				# Solid angle subtended by source on pupil
	A,					# Collecting area of pupil
	eps = 1.0,			# Emissivity of source
	eta = 1.0			# Efficiency of system
	):

	"""
		Given a blackbody with emissivity eps and temperature T, find the total irradiance over an interval 
		[wavelength_min, wavelength_max] incident upon a pupil/emitter system with etendue A * Omega and 
		efficiency eta. 

		I is in units of photons per second (per area A) if eta does not include the QE.
		I is in units of electrons per second (per area A) if eta does include the QE.
	"""

	# Planck function
	B = lambda wavelength_m, T: 2 * constants.h * np.power(constants.c,2) / np.power(wavelength_m,5) * 1 / (np.exp(constants.h * constants.c / (wavelength_m * constants.Boltzmann * T)) - 1)
	
	# Integrand
	integrand = lambda wavelength_m, Omega, eta, A, T, eps: Omega * A * eta * B(wavelength_m, T) * wavelength_m / (constants.h * constants.c) * eps
	
	# Integrating to find the total irradiance incident upon the pupil
	if hasattr(eps, '__call__'):
		# if the emissivity is a function of wavelength_m
		integrand_prime = lambda wavelength_m, Omega, eta, A, T: integrand(wavelength_m, Omega, eta, A, T, eps(wavelength_m))
		I = integrate.quad(integrand_prime, wavelength_min, wavelength_max, args=(Omega, eta, A, T))[0]
	else:	
		# if the emissivity is scalar
		I = integrate.quad(integrand, wavelength_min, wavelength_max, args=(Omega, eta, A, T, eps))[0]

	return I

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
