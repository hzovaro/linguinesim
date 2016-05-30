# Required packages
import scipy.constants as constants
import scipy.integrate as integrate
import numpy as np

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
	B = lambda wavelength, T: 2 * constants.h * np.power(constants.c,2) / np.power(wavelength,5) * 1 / (np.exp(constants.h * constants.c / (wavelength * constants.Boltzmann * T)) - 1)
	
	# Integrand
	integrand = lambda wavelength, Omega, eta, A, T, eps: Omega * A * eta * B(wavelength, T) * wavelength / (constants.h * constants.c) * eps
	
	# Integrating to find the total irradiance incident upon the pupil
	if hasattr(eps, '__call__'):
		# if the emissivity is a function of wavelength
		integrand_prime = lambda wavelength, Omega, eta, A, T: integrand(wavelength, Omega, eta, A, T, eps(wavelength))
		I = integrate.quad(integrand_prime, wavelength_min, wavelength_max, args=(Omega, eta, A, T))[0]
	else:	
		# if the emissivity is scalar
		I = integrate.quad(integrand, wavelength_min, wavelength_max, args=(Omega, eta, A, T, eps))[0]

	return I