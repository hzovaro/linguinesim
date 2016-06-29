############################################################################################
#
# 	File:		sandbox.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#	Edited:		29/06/2016
#
#	Description:
#	A 'sandbox' for testing various features. Play safe kids!
#
#	TO DO:
#	- decide whether ot use classes for telescopes, detectors, etc.
#	- can we pass module names as function args?
#	- tidy up ETC: make a nonspecific function
#
############################################################################################

# A module in Python is like a specialised dictionary that can store functions and which can be 
# variables which can be accessed using the . operator.

from apdsim import *
import apdsim.sysparams.cryo as cryo
import apdsim.sysparams.detector_saphira as detector
import apdsim.sysparams.telescope_anu23m as telescope

# Simulating an image of a galaxy with a Sérsic profile.
# Disc: truncated exponential.
# Bulge: de Vaucoeleurs profile. 

# Inputs:
R_e = 10	# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
i = 45		# inclination (degrees)
alpha = 30	# rotation (degrees)
mu_e = 12	# Surface brightness at the half-light radius (AB mag/arcsec^2)
