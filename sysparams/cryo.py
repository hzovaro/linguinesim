############################################################################################
#
# 	File:		cryoParameters.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Properties of the cryostat that will house the SAPHIRA APD.
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
#	apdsim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with apd-sim.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################################################
from __future__ import division
from linguinesim.sysparams import *

""" Efficiency properties """
Tr_win = 0.98	# Transmission of cryostat window

""" Geometry """
Omega = np.pi 	# Solid angle subtended by the cryostat walls on the detector

""" Thermal properties """
T = 172.372					# Temperature
eps_wall = 1.0			# Wall emissivity
eps_win = 1 - Tr_win	# Window emissivity
print 'TODO: Cryostat temperature needs updating!'