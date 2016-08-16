####################################################################################################
#
# 	File:		etcutils.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Generic utilities for exposure time calcuations.
#
#	Copyright (C) 2016 Anna Zovaro
#
####################################################################################################
#
#	This file is part of lingiune-sim.
#
#	lingiune-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	lingiune-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with lingiune-sim.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################
from __future__ import division
from apdsim import *

def surfaceBrightness2countRate(mu, A_tel, 
	plate_scale_as_px = 1,
	tau = 1,
	qe = 1,
	gain = 1,
	magnitudeSystem = None,
	wavelength_m = None, 
	bandwidth_m = None, 
	band = None
	):
	""" 
		Return the electron count (e/pixel/s) from a source with a given surface brightness OR magnitude (mu) imaged through a system with collecting area A_tel, throughput tau, quantum efficiency qe, internal gain gain and a detector plate scale (plate_scale_as) in arcsec/pixel.

		If mu is given in magnitudes, then the plate scale is irrelevant to the calculation.
	"""
	if (band == None and (wavelength_m == None or bandwidth_m == None)) or (band != None and (wavelength_m != None or bandwidth_m != None)):
		print 'ERROR: you must specify either a band (J, H or K) OR a wavelength_m and bandwidth!'
		return

	if band != None:
		wavelength_m = FILTER_BANDS_M[band][0]
		bandwidth_m = FILTER_BANDS_M[band][1]

	# Getting the magnitude zero points.
	if magnitudeSystem == 'AB':
		zeropoint = AB_MAGNITUDE_ZEROPOINT
	elif magnitudeSystem == 'VEGA':
		if band != None:
			zeropoint = VEGA_MAGNITUDE_ZEROPOINT[band]
		else:
			print 'ERROR: if Vega magnitudes are specified you must also specify the band!'
	else:
		zeropoint = 0

	F = surfaceBrightness2flux(mu=mu, wavelength_m=wavelength_m, zeropoint=zeropoint)
	Sigma_electrons = flux2countRate(F=F, A_tel=A_tel, plate_scale_as_px=plate_scale_as_px, tau=tau, qe=qe, gain=gain, magnitudeSystem=magnitudeSystem, wavelength_m=wavelength_m,bandwidth_m=bandwidth_m, band=band)                                                             
	return Sigma_electrons

####################################################################################################
# This routine is independent of telescope, detector geometry.
def surfaceBrightness2flux(mu, 
	wavelength_m = None,
	zeropoint = AB_MAGNITUDE_ZEROPOINT	# Default: mu is given in AB magnitudes
	):
	"""" 
		Convert a given surface brightness (expressed in magnitudes/arcsec^2) OR magnitude (in magnitudes) to spectral radiance units in both per unit frequency and per unit wavelength_m.
		The flux values are returned in both CGS and SI units.

		If a wavelength_m is not specified then only the spectral radiance per unit 
		frequency is returned.
	"""
	F_nu_cgs = np.power(10, - (zeropoint + mu) / 2.5)						# ergs/s/cm^2/arcsec^2/Hz
	return Fnucgs2flux(F_nu_cgs, wavelength_m)


####################################################################################################
def Fnucgs2flux(F_nu_cgs,
	wavelength_m = None
	):
	"""
		Convert a given flux frequency density specified in cgs units (ergs/s/Hz/cm^-2) to cgs flux wavelength density units (ergs/s/Angstrom/cm^-2). Also returns these two flux values in SI units.
	"""

	F_nu_si = F_nu_cgs * 1e-7 * 1e4 # W/m^2/Hz
	if wavelength_m != None:
		F_lambda_cgs = F_nu_cgs * constants.c / np.power(wavelength_m, 2) * 1e-10	# ergs/s/cm^2/angstrom
		F_lambda_si = F_lambda_cgs * 1e-7 * 1e4	* 1e10								# W/m^2/m
	else:
		F_lambda_cgs = None
		F_lambda_si = None

	F = {
		'F_nu_cgs' 		: F_nu_cgs,
		'F_lambda_cgs' 	: F_lambda_cgs,
		'F_nu_si' 		: F_nu_si,
		'F_lambda_si' 	: F_lambda_si,
	}

	return F

####################################################################################################
def flux2photonRate(F, wavelength_m, bandwidth_m):
	"""
		Convert a given flux from a source (in a dictionary format output by surfaceBrightness2flux) into photons/s/m^2/arcsec^2 (or in photons/s/m^2 depending on the units of F) given a central wavelength_m and bandwidth of a filter.
	"""
	E_photon = constants.h * constants.c / wavelength_m			# J
	Sigma_photons = F['F_lambda_si'] * bandwidth_m / E_photon	# W/m^2/arcsec^2/m * m/J = photons/s/m^2/arcsec^2
	return Sigma_photons

####################################################################################################
def photonRate2countRate(Sigma_photons, A_tel, 
	tau, 
	qe, 
	gain,
	plate_scale_as_px = 1
	):
	"""
		Convert a flux given in units of photons/s/m^2/arcsec^2 (or in units of photons/s/m^2 depending on the units of Sigma_photons) into detector counts given a telescope collecting area (A_tel), detector plate scale (arcsec/pixel), throughput (tau), quantum efficiency (qe) and gain (gain). 

		If the input is not specified per square arcsec, then the plate_scale_as_px should be left blank.
	"""
	Sigma_electrons = Sigma_photons * A_tel * plate_scale_as_px * plate_scale_as_px * tau * qe * gain # photons/s/px
	return Sigma_electrons

####################################################################################################
def flux2countRate(F, A_tel, 
	plate_scale_as_px = 1, 
	tau = 1,
	qe = 1,
	gain = 1,
	magnitudeSystem = None,
	wavelength_m = None, 
	bandwidth_m = None, 
	band = None
	):
	""" 
		Return the electron count (e/pixel/s) from a source with a given 
		imaged through a system with collecting area A_tel, throughput tau, quantum efficiency qe, 
		internal gain gain and a detector plate scale (plate_scale_as) in arcsec/pixel 
	"""

	if band != None:
		wavelength_m = FILTER_BANDS_M[band][0]
		bandwidth_m = FILTER_BANDS_M[band][1]
	elif (wavelength_m == None or bandwidth_m == None):
		print 'ERROR: you must specify either a band (J, H or K) OR a wavelength_m and bandwidth!'
		return

	Sigma_photons = flux2photonRate(F=F, wavelength_m=wavelength_m, bandwidth_m=bandwidth_m)
	Sigma_electrons = photonRate2countRate(Sigma_photons=Sigma_photons, A_tel=A_tel, plate_scale_as_px=plate_scale_as_px, tau=tau, qe=qe, gain=gain)
	# pdb.set_trace()
	return Sigma_electrons

####################################################################################################
def expectedCount2count(arg,
	detectorSaturation = np.inf, 
	t_exp = None):
	"""
		Convert an expected photon count (in photons) OR expected photon count rate (in photons/s) to a 'truth' count assuming a Poisson distribution.
		If t_exp is not specified, it is assumed that the input argument is given in units of photons. 
		If t_exp is specified, it is assumed that the input argument is given in units of photons/s, in which case the expected count is arg * t_exp.
		The saturation level of the detector may also be specified if desired. By default this value is infinity.
	"""
	if t_exp == None:
		expectedCount = arg
	else:
		expectedCount = arg * t_exp
	return np.random.poisson(lam=np.minimum(np.maximum(expectedCount, 0), detectorSaturation), size=expectedCount.shape)