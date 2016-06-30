# Taken from http://www.astro.ucla.edu/%7Ewright/CosmoCalc.html 

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt

###################################################################################
def distances(z=3.0, verbose=False):

	H0 = 70                        # Hubble constant
	WM = 0.3                       # Omega(matter)
	WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

	# initialize constants
	WR = 0.        # Omega(radiation)
	WK = 0.        # Omega curvaturve = 1-Omega(total)
	c = 299792.458 # velocity of light in km/sec
	Tyr = 977.8    # coefficent for converting 1/H into Gyr
	DTT = 0.5      # time from z to now in units of 1/H0
	DTT_Gyr = 0.0  # value of DTT in Gyr
	age = 0.5      # age of Universe in units of 1/H0
	age_Gyr = 0.0  # value of age in Gyr
	zage = 0.1     # age of Universe at redshift z in units of 1/H0
	zage_Gyr = 0.0 # value of zage in Gyr
	DCMR = 0.0     # comoving radial distance in units of c/H0
	DCMR_Mpc = 0.0 
	DCMR_Gyr = 0.0
	DA = 0.0       # angular size distance
	DA_Mpc = 0.0
	DA_Gyr = 0.0
	kpc_DA = 0.0
	DL = 0.0       # luminosity distance
	DL_Mpc = 0.0
	DL_Gyr = 0.0   # DL in units of billions of light years
	V_Gpc = 0.0
	a = 1.0        # 1/(1+z), the scale factor of the Universe
	az = 0.5       # 1/(1+z(object))

	h = H0/100.
	WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
	WK = 1-WM-WR-WV
	az = 1.0/(1+1.0*z)
	age = 0.
	n=1000         # number of points in integrals

	for i in range(n):
		a = az*(i+0.5)/n
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
		age = age + 1./adot

	zage = az*age/n
	zage_Gyr = (Tyr/H0)*zage
	DTT = 0.0
	DCMR = 0.0

	# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
	for i in range(n):
		a = az+(1-az)*(i+0.5)/n
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
		DTT = DTT + 1./adot
		DCMR = DCMR + 1./(a*adot)

	DTT = (1.-az)*DTT/n
	DCMR = (1.-az)*DCMR/n
	age = DTT+zage
	age_Gyr = age*(Tyr/H0)
	DTT_Gyr = (Tyr/H0)*DTT
	DCMR_Gyr = (Tyr/H0)*DCMR
	DCMR_Mpc = (c/H0)*DCMR

	# tangential comoving distance

	ratio = 1.00
	x = sqrt(abs(WK))*DCMR
	if x > 0.1:
		if WK > 0:
		  ratio =  0.5*(exp(x)-exp(-x))/x 
		else:
		  ratio = sin(x)/x
	else:
		y = x*x
		if WK < 0: y = -y
		ratio = 1. + y/6. + y*y/120.
	DCMT = ratio*DCMR
	DA = az*DCMT
	DA_Mpc = (c/H0)*DA
	kpc_DA = DA_Mpc/206.264806
	DA_Gyr = (Tyr/H0)*DA
	DL = DA/(az*az)
	DL_Mpc = (c/H0)*DL
	DL_Gyr = (Tyr/H0)*DL

	# comoving volume computation

	ratio = 1.00
	x = sqrt(abs(WK))*DCMR
	if x > 0.1:
		if WK > 0:
		  ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
		else:
		  ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
	else:
		y = x*x
		if WK < 0: y = -y
		ratio = 1. + y/5. + (2./105.)*y*y
	VCM = ratio*DCMR*DCMR*DCMR/3.
	V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM

	if verbose:
		print 'For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ',
		print '%1.2f' % WV + ', z = ' + '%1.3f' % z
		print 'It is now ' + '%1.1f' % age_Gyr + ' Gyr since the Big Bang.'
		print 'The age at redshift z was ' + '%1.1f' % zage_Gyr + ' Gyr.'
		print 'The light travel time was ' + '%1.1f' % DTT_Gyr + ' Gyr.'
		print 'The comoving radial distance, which goes into Hubbles law, is',
		print '%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % DCMR_Gyr + ' Gly.'
		print 'The comoving volume within redshift z is ' + '%1.1f' % V_Gpc + ' Gpc^3.'
		print 'The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or',
		print '%1.1f' % DA_Gyr + ' Gly.'
		print 'This gives a scale of ' + '%.2f' % kpc_DA + ' kpc/".'
		print 'The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc or ' + '%1.1f' % DL_Gyr + ' Gly.'
		print 'The distance modulus, m-M, is '+'%1.2f' % (5*log10(DL_Mpc*1e6)-5)

	return {
		'D_now_Mpc' : DCMR_Mpc,
		'D_A_Mpc'	: DA_Mpc,
		'D_L_Mpc'	: DL_Mpc
	}

###################################################################################
# def thetaVsRedshift():
" Calculate the angular size subtended on the size by 1 kpc as a function of redshift. "
# Diffraction-limited resolution at the ANU 2.3 m
res_ANU_J = 1.25e-6 / 2.3 * 180 / np.pi * 3600
res_ANU_H = 1.635e-6 / 2.3 * 180 / np.pi * 3600
res_ANU_K = 2.2e-6 / 2.3 * 180 / np.pi * 3600
# For SAMI
res_SAMI = 1.6

z = np.linspace(0.004,0.095,500)
D_A = np.zeros(z.shape)
D_L = np.zeros(z.shape)
D_now = np.zeros(z.shape) 
# theta_arc = np.zeros(z.shape)
theta_D_A = np.zeros(z.shape)
theta_D_L = np.zeros(z.shape)
theta_D_now = np.zeros(z.shape)

for k in range(z.size):
	dist = distances(z[k])
	D_A[k] = dist['D_A_Mpc']
	D_L[k] = dist['D_L_Mpc']
	D_now[k] = dist['D_now_Mpc']
	theta_D_A[k] = 2 * np.arctan(1e-3 / 2 / D_A[k]) * 180 / np.pi * 3600

# Plotting
plt.close('all')
plt.figure()
plt.subplot(1,2,1)
plt.semilogy(z, theta_D_A, 'g', label=r'$\theta$')
# SAMI resolution as reference
plt.semilogy(z, np.ones(z.shape) * res_ANU_J, color='darkorange', label=r'2.3 m, $J$-band')
plt.semilogy(z, np.ones(z.shape) * res_ANU_H, color='orangered', label=r'2.3 m, $H$-band')
plt.semilogy(z, np.ones(z.shape) * res_ANU_K, color='darkred', label=r'2.3 m, $K$-band')
plt.semilogy(z, np.ones(z.shape) * 1.6, 'k--', label='SAMI (diffraction-limited)')
plt.semilogy(z, np.ones(z.shape) * 2.1, linestyle='dashed', color='darkgrey', label='SAMI (median seeing)')
plt.semilogy(z, np.ones(z.shape) * 1., linestyle='dashed', color='blue', label='Ground-based imaging data resolution')
plt.xlim([min(z), max(z)])
plt.legend()
# plt.xlabel(r'$D_A(z)$ (Mpc)')
plt.xlabel(r'Redshift $z$')
# plt.ylabel(r'Angular size $\log_{10}(\theta)$ (arcsec)')
plt.ylabel(r'Angular size $\theta$ (arcsec)')
plt.title(r'Angular size subtended on sky by 1 kpc')

" Calculate the distance correspoding to a given on-sky angle as a function of redshift "
# SAMI
l_groundbased = 2 * np.tan(1. / 3600 * np.pi / 180 / 2) * D_A * 1e3
l_SAMI = 2 * np.tan(1.6 / 3600 * np.pi / 180 / 2) * D_A * 1e3
l_SAMI_median = 2 * np.tan(2.1 / 3600 * np.pi / 180 / 2) * D_A * 1e3
l_ANU_J = 2 * np.tan(res_ANU_J / 3600 * np.pi / 180 / 2) * D_A * 1e3
l_ANU_H = 2 * np.tan(res_ANU_H / 3600 * np.pi / 180 / 2) * D_A * 1e3
l_ANU_K = 2 * np.tan(res_ANU_K / 3600 * np.pi / 180 / 2) * D_A * 1e3

plt.subplot(1,2,2)
plt.plot(z, l_ANU_J, color='darkorange', label=r'2.3 m, $J$-band')
plt.plot(z, l_ANU_H, color='orangered', label=r'2.3 m, $H$-band')
plt.plot(z, l_ANU_K, color='darkred', label=r'2.3 m, $K$-band')
plt.plot(z, l_SAMI, 'k--', label='SAMI (diffraction-limited)')
plt.plot(z, l_SAMI_median, linestyle='dashed', color='darkgrey', label='SAMI (median seeing)')
plt.plot(z, l_groundbased, linestyle='dashed', color='blue', label='Ground-based imaging data resolution')
plt.xlim([min(z), max(z)])
plt.legend(loc='upper left')
plt.xlabel(r'Redshift $z$')
plt.ylabel(r'$l$ (kpc)')
plt.title(r'Distance corresponding to diffraction-limited resolution')
plt.show()

" Cosmological distances "
# plt.figure()
# plt.plot(z, D_A, 'r', label=r'$D_A(z)$')
# plt.plot(z, D_L, 'b', label=r'$D_L(z)$')
# plt.plot(z, D_now, 'g', label=r'$D_{now}(z)$')
# plt.xlim([min(z), max(z)])
# plt.legend(loc='lower right')
# plt.xlabel(r'Redshift $z$')
# plt.ylabel(r'Distance (Mpc)')
# plt.title(r'Cosmological distances')
# plt.show()
