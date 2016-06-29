############################################################################################
#
# 	File:		h4rgParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		17/06/2016
#
#	Description:
#	Properties of the HAWAII H4RG detector.
#
############################################################################################
from sysparams import *

""" Efficiency properties """
gain_av = 1					# avalanche gain (electrons/photelectron)
qe = 1.0					# quantum efficiency (worst-case)
adu_gain = 1 / 2.9			# electrons/ADU 

""" Geometry """
width_px = 4096				# width (pixels)
height_px = 4096			# height (pixels)
n_px_tot = width_px * height_px	# number of pixels (total)
l_px = 10e-6				# pixel width (m)
# l_px = 15e-6				# pixel width (m)
A_px = l_px * l_px			# pixel area (m^2)

""" Noise properties """
read_noise = 16				# sqrt(e/pixel) rms (CDS)
dark_current = 0.05			# e/second/pixel. Note that the dark current will be smaller for 10um pixels. Source: http://www.gmto.org/wp-content/uploads/TELEDYNE%20Detector%20Update%20(for%20ELTs)%20-%2023%20Oct%202015.pdf
wavelength_cutoff = 1.8e-6	# cutoff wavelength (m) - 'K-blind'
wavelength_cutoff_h = 3e-6	# worst-case cutoff wavelength (m)