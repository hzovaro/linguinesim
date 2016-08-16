#########################################################################################################
#
# 	File:		linguini_unittest.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Unit testing of the etc, etcutils, galsim, imutils, lisim and obssim modules.
#
#	Copyright (C) 2016 Anna Zovaro
#
#########################################################################################################
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
#########################################################################################################
from __future__ import division
import pprint
pp = pprint.PrettyPrinter(indent=4)
from apdsim import *

def etcutils_test():
	print "TEST: surfaceBrightness2flux()"
	# Input parameters:
	mu = 15
	wavelength = 900e-9
	
	# Outputs (produced):
	F = surfaceBrightness2flux(mu, wavelength)
	pp.pprint(F)

	print "\nTEST: flux2photonRate()"
	# Input parameters:
	bandwidth = 10e-9
	# Outputs (produced):
	Sigma_photons = flux2photonRate(F, wavelength, bandwidth)
	print "F_phot:", Sigma_photons

	print "\nTEST: photonRate2countRate()"
	# Input parameters:
	Sigma_photons = 600 * 1e4
	A_tel = np.pi * 4 * 4
	plate_scale = 0.45 
	tau = 0.98
	qe = 0.91
	gain = 60
	# Outputs (produced):
	Sigma_electrons = photonRate2countRate(Sigma_photons, A_tel, plate_scale, tau, qe, gain)
	print "F_electrons:", Sigma_electrons

	print "\nTEST: flux2countRate()"
	# Input parameters: 
	# Outputs (produced):
	Sigma_electrons = flux2countRate(F, A_tel, plate_scale, tau, qe, gain, magnitudeSystem='AB', wavelength_m=wavelength, bandwidth_m=bandwidth)
	print "F_electrons:", Sigma_electrons
	return


