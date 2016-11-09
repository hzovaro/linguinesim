################################################################################
#
# 	File:		galaxy.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A class for a galaxy.
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

################################################################################
class Galaxy(object):

	def __init__(self, name, R_e_as, mu_e, sersic_idx,
		axis_ratio=1,
		PA_deg=0,
		z=None,
		masses_solar=None,
		coords=None,
		gama_id=None):

		self.name = name
		self.R_e_as = R_e_as
		self.mu_e = mu_e
		self.sersic_idx = sersic_idx
		self.axis_ratio = axis_ratio
		self.PA_deg = PA_deg
		self.z = z
		self.masses_solar = masses_solar

		self.coords = coords
		self.gama_id = gama_id