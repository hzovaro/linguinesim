############################################################################################
#
# 	File:		sandbox.py
#	Author:		Anna Zovaro
#	Email:		anna.zovaro@anu.edu.au
#
#	Description:
#	A 'sandbox' for testing various features. Play safe kids!
#
#	TO DO:
#	- Check consistency of np.log usage!
#
############################################################################################
from apdsim import *

def sersic(n, R_e, 
	n_points=100, 
	R_max = None,
	mu_e = None,
	I_e = None
	):
	" Returns a one-dimensional Sersic profile in units of surface brightness with index n, half-light radius R_e and mu(R=R_e) = mu_e. "
	
	if R_max == None:
		R_max = 10 * R_e
	R = np.linspace(-R_max, +R_max, n_points)
	
	# Calculating the constraining parameter.
	if mu_e == None and I_e == None:
		print 'ERROR: at least one of mu_e or I_e must be specified!'
		return
	if mu_e == None:
		mu_e = - 2.5 * np.log10(I_e)
	elif I_e == None:
		I_e = np.power(10, - mu_e / 2.5)
	else:
		print 'WARNING: you have specified both mu_e and I_e, so they might not be consistent! Proceed with caution...'

	# Calculating b_n given the Sersic index n.
	if n > 0.5 and n < 8:
		b_n = 1.9992 * n - 0.3271
	elif n >= 8:
		b_n = 2 * n - 1 / 3
	else:
		print 'ERROR: haven\'t implemented n < 0.5 yet! Check again later...'
		return

	mu = mu_e + 2.5 / np.log(10) * b_n * (np.power(R/R_e, n) - 1)
	# NOTE: I is in units of ergs/s/cm^2/arcsec^2/Hz? Check with Rob...
	I = I_e * np.exp(- b_n * (np.power(R/R_e, 1/n)))

	return R, mu, I


# Simulating an image of a galaxy with a Sersic profile.
# Disc: truncated exponential.
# Bulge: de Vaucoeleurs profile. 

" Inputs "
# Galaxy:
R_e = 10	# Half-light or effective radius (kpc) - that is, the radius enclosing half the total light from the galaxy
mu_e = 12	# Surface brightness at the half-light radius (AB mag/arcsec^2)
n = 2
z = 0.001	# Redshift
i = 45		# inclination (degrees)
alpha = 30	# rotation (degrees)

# So the galaxy has a certain redshift and a certain half-light radius
# We want our process to be: Light profile in kpc --> Light profile in arcsec --> Light profile in pixels.
# So make the Sersic profile in units of linear distance (say pc or kpc), then use the redshift to convert 1 kpc into arcseconds acros the sky (use cosmo_calc.py for this), then rescale the Sersic profile to arcseconds, then rescale to pixels on the detector.
R, mu, I = sersic(n=n, R_e=10)

plt.figure()
plt.subplot(1,2,1)
plt.plot(R, mu, 'r')
plt.title(r'Surface brightness vs. radius, $n = %d$' % n)
plt.xlabel(r'Radius $R$ (kpc)')
plt.ylabel(r'Surface brightness $\mu(R)$ (magnitudes/arcsec$^2$)')
plt.subplot(1,2,2)
plt.plot(R, mu, 'b')
plt.title(r'Intensity vs. radius, $n = %d$' % n)
plt.xlabel(r'Radius $R$ (kpc)')
plt.ylabel(r'Intensity $I(R)$')
