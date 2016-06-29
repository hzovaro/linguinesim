############################################################################################
#
# 	File:		anu23mParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		27/06/2016
#
#	Description:
#	Properties of the ANU 2.3 m telescope at Siding Spring Observatory.
#	The telescope is assumed to be used in the nasmyth configuration.
#
###########################################################################################
from scipy import constants
import numpy as np

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
A_collecting = A_M1_reflective

""" f ratio & plate scale """
efl_mm = 6010. / 700. * 4840.	# effective focal length
f_ratio = efl_mm / 1e3 / 2 / r_M1_outer
plate_scale_as_mm = 206256. / efl_mm
plate_scale_as_m = plate_scale_as_mm * 1e3 

""" Imaging filter bands """
filter_bands_um = {
	# [centre wavelength, width, min, max]
	# 'Y' : [1.030, 0.100],	
	'J' : [1.250, 0.160, 1.170, 1.330],
	'H' : [1.635, 0.290, 1.490, 1.780],
	'K' : [2.200, 0.340, 2.030, 2.370]
}

filter_bands_m = {
	# [centre wavelength, width, min, max]
	'J' : [1.250e-6, 0.160e-6, 1.170e-6, 1.330e-6],
	'H' : [1.635e-6, 0.290e-6, 1.490e-6, 1.780e-6],
	'K' : [2.200e-6, 0.340e-6, 2.030e-6, 2.370e-6]
}

""" Sky brightness (magnitudes per square arcsec) """
# Source: http://www.mso.anu.edu.au/pfrancis/reference/reference/node4.html
# Note these are 'full moon' values -- typical values will be better!
sky_brightness = {
	'J' : 15,
	'H' : 13.7,
	'K' : 12.5
}