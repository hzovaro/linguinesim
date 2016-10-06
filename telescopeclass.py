############################################################################################
#
# 	File:		telescopeclass.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A class for a telescope.
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
from linguineglobals import *

####################################################################################################
class Telescope(object):

	def __init__(self,
		efl_m,
		T=273.15
		):

		self.T = T
		self.tau = 1.0	# initialise the throughput to 1 as it is determined by the telescope's mirrors.
		self.efl_m = efl_m
		self.efl_mm = efl_m * 1e3
		
		self.f_ratio = None

		self.plate_scale_as_mm = 206265 / self.efl_mm
		self.plate_scale_as_m = 206265 / self.efl_mm * 1e3
		self.plate_scale_rad_mm = self.plate_scale_as_mm / 3600 * np.pi / 180
		self.plate_scale_rad_m = self.plate_scale_as_m / 3600 * np.pi / 180

		# Initialise the mirror list
		self.mirrors = []
		self._N_mirrors = 0

		# Spider
		self.hasSpider = False

	def addMirror(self,
		R_outer_m,
		R_inner_m=0.0,
		reflectivity=AL_REFLECTIVITY,
		eps_eff=None
		):
		
		newmirror = MirrorClass(
			R_outer_m=R_outer_m,
			R_inner_m=R_inner_m,
			reflectivity=reflectivity,
			eps_eff=eps_eff)

		# Append to the telescope's list of mirrors.	
		self.mirrors.append(newmirror)
		self._N_mirrors = len(self.mirrors)

		# We set the f ratio using the diameter of M1.
		if self._N_mirrors == 1:
			self.A_collecting_m2 = self.mirrors[0].A_reflective_m2
			self.f_ratio = self.efl_m / self.mirrors[0].D_outer_m

		# Update the throughput.
		self.tau *= newmirror.reflectivity

	def addSpider(self,
		A_spider_m2,
		eps_spider=AL_EMISSIVITY):

		self.hasSpider = True
		# The spider will be eps_spider emissive at telescope temp, and (1 - eps_spider) * eps_sky emissive at the sky temp.
		self.A_spider_m2 = A_spider_m2
		self.eps_spider = eps_spider
		self.eps_spider_eff = A_spider_m2 / self.A_collecting_m2 * self.eps_spider

####################################################################################################
class MirrorClass(object):

	def __init__(self,
		R_outer_m,
		R_inner_m,
		reflectivity=AL_REFLECTIVITY,
		eps_eff=None
		):
		"""
			Create a mirror object. Mirrors are defined by an outer and inner radius, a reflectivity and an effective emissivity.

			The effective emissivity can be specified manually, or by default. If unspecified, the effective emissivity is the mean emissivity of the mirror weighted by the area of the central hole and by that of the reflective part. 
		"""

		# Geometry
		self.R_outer_m = R_outer_m
		self.R_inner_m = R_inner_m
		self.D_outer_m = self.R_outer_m * 2
		self.D_inner_m = self.R_inner_m * 2
		self.A_hole_m2 = np.pi * self.R_inner_m**2
		self.A_outer_m2 = np.pi * self.R_outer_m**2	# use this for ETC calculations
		self.A_reflective_m2 = self.A_outer_m2 - self.A_hole_m2

		# Thermal characteristics
		self.reflectivity = reflectivity
		self.eps_reflective = 1 - self.reflectivity
		self.eps_hole = 1.0

		# If the effective emissivity is not specified, then we calculate it assuming 
		if not eps_eff:
			self.eps_eff = (self.A_hole_m2 * self.eps_hole + self.A_reflective_m2 * self.eps_reflective) / self.A_outer_m2
		else:
			self.eps_eff = eps_eff
	