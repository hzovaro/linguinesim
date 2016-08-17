############################################################################################
#
# 	File:		sky_sso.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	Sky brightness and temperature properties recorded at Siding Spring Observatory.
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
from linguinesim.sysparams import *

""" Sky brightness (magnitudes per square arcsec) """
# Source: http://www.mso.anu.edu.au/pfrancis/reference/reference/node4.html
# Note these are 'full moon' values -- typical values will be better!
magnitude_system = 'AB'
brightness = {
	'J' : 15,
	'H' : 13.7,
	'K' : 12.5
}
T = 273