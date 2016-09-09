############################################################################################
#
# 	File:		detectorclass.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A class for a cryostat.
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
from __future__ import division, print_function
import numpy as np

####################################################################################################
class Detector(object):

	def __init__(self,
		height_px,
		width_px,
		l_px_m,
		wavelength_cutoff=np.inf,
		wavelength_cutoff_h=np.inf,
		RN=0,
		dark_current=0,
		gain=1,
		saturation=np.inf,
		adu_gain=1,
		qe=1):

		# Detector geometry.
		self.height_px = height_px
		self.width_px = width_px
		self.size_px = (height_px, width_px)
		self.l_px_m = l_px_m
		self.A_px_m2 = l_px_m**2

		# Electronic properties.		
		self.gain = gain 	# Internal gain
		self.qe = qe
		self.adu_gain = adu_gain
		self.saturation = saturation
		self.RN = RN

		# Optical properties.
		self.wavelength_cutoff = wavelength_cutoff
		self.wavelength_cutoff_h = wavelength_cutoff_h

		# Noise properties.
		self.RN = RN
		self.dark_current = dark_current