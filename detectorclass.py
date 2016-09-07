from __future__ import division, print_function
from linguinesim.apdsim import *

class DetectorClass(object):

	def __init__(self,
		height_px,
		width_px,
		l_px_m,
		wavelength_cutoff=np.inf,
		wavelength_cutoff_h=np.inf,
		RN=0,
		dark_current=0,
		gain=1,
		saturation=np.inf,
		adu_gain=1,
		qe=1):

		# Detector geometry.
		self.height_px = height_px
		self.width_px = width_px
		self.size_px = (height_px, width_px)
		self.l_px_m = l_px_m
		self.A_px_m2 = l_px_m**2

		# Electronic properties.		
		self.gain = gain 	# Internal gain
		self.qe = qe
		self.adu_gain = adu_gain
		self.saturation = saturation
		self.RN = RN

		# Optical properties.
		self.wavelength_cutoff = wavelength_cutoff
		self.wavelength_cutoff_h = wavelength_cutoff_h

		# Noise properties.
		self.RN = RN
		self.dark_current = dark_current