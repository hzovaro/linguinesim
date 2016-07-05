############################################################################################
#
# 	File:		anu23mParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Properties of the ANU 2.3 m telescope at Siding Spring Observatory.
#	The telescope is assumed to be used in the nasmyth configuration.
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

""" Mirror geometry & thermal properties """
T = 273 + 10		# telescope temperature (K)

""" Primary mirror """
r_M1_outer = 2.337 / 2	# primary mirror radius (m) (from http://rsaa.anu.edu.au/observatories/telescopes/anu-23m-telescope)
r_M1_inner = 0.229 / 2	# central hole radius (m)

A_M1_hole = np.power(r_M1_inner,2) * np.pi 		# Primary (central hole) (m^2)
A_M1_total = np.power(r_M1_outer,2) * np.pi 	# Primary (total area) (m^2)
A_M1_reflective = A_M1_total - A_M1_hole		# Primary (reflective) (m^2)

eps_M1_reflective = 0.091		# emissivity of reflective part
eps_M1_hole = 1.0				# emissivity of central hole
eps_M1_eff = A_M1_hole / A_M1_total * eps_M1_hole + A_M1_reflective / A_M1_total * eps_M1_reflective

""" Secondary mirror (Nasmyth) """
# Note: using the worst-case value corresponding to the baffle radius
# ACTUAL radii
r_M2_outer = 0.355 / 2						# ACTUAL radius of M2
r_M2_baffle = 0.500 / 2						# ACTUAL radius of baffle
r_M2_oversize = r_M2_baffle - r_M2_outer 	# ACTUAL annular radius of oversized-ness 

# EFFECTIVE radii (scaled up to M1 size)
r_M2_oversize_eff = r_M2_oversize / r_M2_baffle * r_M1_outer			# EFFECTIVE oversized-ness radius of M2 (scaled up to M1)
r_M2_outer_eff = r_M1_outer + r_M2_oversize_eff						# EFFECTIVE radius of M2 (scaled up to M1)

A_M2_total_eff = np.power(r_M2_outer_eff, 2) * np.pi
A_M2_reflective_eff = A_M1_total
A_M2_oversize_eff = A_M2_total_eff - A_M1_total		# EFFECTIVE oversized area of M2 (scaled up to M1)

eps_M2_reflective = 0.091		# emissivity of reflective part
eps_baffle = 1.0	# emissivity of oversized portion
eps_M2_eff = A_M2_oversize_eff / A_M2_total_eff * eps_baffle + A_M2_reflective_eff / A_M2_total_eff * eps_M2_reflective

""" Tertiary """
eps_M3_reflective = 0.091

""" Spider """
A_spider = 4 * 0.012 * (r_M1_outer - r_M2_outer)	# Total area of spider (m^2)
eps_spider = 0.091	# emissivity of spider
eps_spider_eff = A_spider / A_M1_total * eps_spider 

""" Telescope throughput """
tau = 0.909 * 0.909 * 0.909	# telescope transmission	

""" More readable parameters """
r_M1 = r_M1_outer
r_M2 = r_M2_outer
D_M1 = 2 * r_M1
A_collecting = A_M1_reflective

""" f ratio & plate scale """
efl_mm = 6010. / 700. * 4840.	# effective focal length
efl_m = efl_mm / 1e3
f_ratio = efl_m / D_M1
plate_scale_as_mm = 206256. / efl_mm
plate_scale_as_m = plate_scale_as_mm * 1e3 

""" Sky brightness (magnitudes per square arcsec) """
# Source: http://www.mso.anu.edu.au/pfrancis/reference/reference/node4.html
# Note these are 'full moon' values -- typical values will be better!
sky_brightness_magnitude_system = 'AB'
sky_brightness = {
	'J' : 15,
	'H' : 13.7,
	'K' : 12.5
}
T_sky = 273