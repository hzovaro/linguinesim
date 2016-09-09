############################################################################################
#
# 	File:		cryostatclass.py
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

####################################################################################################
class Cryostat(object):

	def __init__(self,
		T, 
		Tr_win,
		Omega, 
		eps_wall
		):

		self.T = T
		self.Tr_win = Tr_win
		self.Omega = Omega
		self.eps_wall = eps_wall
		self.eps_win = 1 - self.Tr_win