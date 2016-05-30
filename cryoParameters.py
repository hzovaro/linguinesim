############################################################################################
#
# 	File:		cryoParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		26/05/2016
#
#	Description:
#	Properties of the cryostat that will house the SAPHIRA APD.
#
############################################################################################
import numpy as np

""" Efficiency properties """
Tr_win = 0.98	# Transmission of cryostat window

""" Geometry """
Omega = np.pi 	# Solid angle subtended by the cryostat walls on the detector

""" Thermal properties """
T = 172.372					# Temperature
eps_wall = 1.0			# Wall emissivity
eps_win = 1 - Tr_win	# Window emissivity