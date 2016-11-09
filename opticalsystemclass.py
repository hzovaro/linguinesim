################################################################################
#
# 	File:		opticalsystemclass.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A class for an optical system incuding a detector, telescope and (optionally) a cryostat and sky.
#
#	Copyright (C) 2016 Anna Zovaro
#
################################################################################
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
################################################################################
from __future__ import division, print_function
import numpy as np

################################################################################
class OpticalSystem(object):

	def __init__(self, telescope, detector, sky,
		plate_scale_as_px=None,
		cryostat=None
		):
		self.telescope = telescope
		self.detector = detector
		self.cryostat = cryostat
		self.sky = sky

		# Plate scales
		if not plate_scale_as_px:
			self.plate_scale_as_px = self.detector.l_px_m * telescope.plate_scale_as_m
			self.plate_scale_rad_px = self.detector.l_px_m * telescope.plate_scale_rad_m
		else:
			self.plate_scale_as_px = plate_scale_as_px
			self.plate_scale_rad_px = np.deg2rad(self.plate_scale_as_px / 3600)
		
		# Detector field of view (FoV)
		self.FoV_height_as = self.detector.height_px * self.plate_scale_as_px
		self.FoV_width_as = self.detector.width_px * self.plate_scale_as_px
		
		self.FoV_height_rad = self.detector.height_px * self.plate_scale_rad_px
		self.FoV_width_rad = self.detector.width_px * self.plate_scale_rad_px
		
		self.FoV_diag_as = np.sqrt(self.FoV_height_as**2 + self.FoV_width_as**2)
		self.FoV_diag_rad = np.sqrt(self.FoV_height_rad**2 + self.FoV_width_rad**2)

		# Pixel FoV
		self.omega_px_as2 = self.plate_scale_as_px**2
		self.omega_px_sr = self.plate_scale_rad_px**2

		# Etendue
		self.etendue = self.omega_px_sr * self.telescope.A_collecting_m2




	