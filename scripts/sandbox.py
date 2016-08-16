############################################################################################
#
# 	File:		sandbox.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A 'sandbox' for testing various features. Play safe kids!
#
#	Copyright (C) 2016 Anna Zovaro
#
############################################################################################
#
#	This file is part of lingiune-sim.
#
#	lingiune-sim is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	lingiune-sim is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with lingiune-sim.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################
from __future__ import division
from apdsim import *
from astropy.modeling import models, fitting
plt.close('all')

# Testing the getDiffractionLimitedImage() function...
pdb.set_trace()
fname = '/Users/azovaro/Documents/Modules/scripts/satellite_images/hst_original.fit'
fitsFileData = apdsim.getRawImages(fname)
im_truth = apdsim.rotateAndCrop(fitsFileData[0])
im_truth = im_truth[0]
width_px = im_truth.shape[0]
height_px = im_truth.shape[1]
hst_length_px = 0.5 * width_px
# Figuring out the plate scale in arcsec...
hst_alt_m = 559e3	# orbital altitude (m)
hst_length_m = 13.2	# HST length (m)
hst_length_as = np.rad2deg(hst_length_m / hst_alt_m) * 3600
hst_as_px = hst_length_as / hst_length_px

# Resize it first so it isn't too big
im_truth = resizeImagesToDetector(images_raw = im_truth, source_plate_scale_as = hst_as_px, dest_plate_scale_as = hst_as_px * 4)
width_px = im_truth.shape[0]
height_px = im_truth.shape[1]
