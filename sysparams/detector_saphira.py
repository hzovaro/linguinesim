############################################################################################
#
# 	File:		apdParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Properties of the SELEX avalanche photodiode detector (APD) with the SAPHIRA ROIC.
#	All properties sourced from Atksinson et al. (2014) - "Observatory Deployment and 
#	Characterization of SAPHIRA HgCdTe APD Arrays"
#
#	Copyright (C) 2016 Anna Zovaro
#
#########################################################################################################
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
#########################################################################################################
from __future__ import division
from sysparams import *

""" Efficiency properties """
gain_av = 52				# avalanche gain (electrons/photelectron)
qe = 0.9					# quantum efficiency (worst-case)
adu_gain = 1 / 2.9			# electrons/ADU 

""" Geometry """
width_px = 320				# width (pixels)
height_px = 256				# height (pixels)
n_px_tot = width_px * height_px	# number of pixels (total)
l_px_m = 24e-6				# pixel width (m)
A_px_m2 = l_px_m * l_px_m	# pixel area (m^2)

""" Noise properties """
read_noise = 9				# sqrt(e/pixel) rms 
# dark_current = (264 + 62)	# e/second/pixel; worst-case non-gain-corrected dark current (see Atkinson et al. 2014)
dark_current = 0.03
wavelength_cutoff = 2.5e-6	# cutoff wavelength (m)
wavelength_cutoff_h = 3e-6	# worst-case cutoff wavelength (m)