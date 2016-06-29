# Required packages
import scipy.constants as constants
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

T = np.array([0., 20.]) + 273.
wavelengths = np.linspace(0.5e-6, 2.5e-6, 1000)
B = np.zeros((len(T), len(wavelengths)))
B_fun = lambda wavelength, T: 2 * constants.h * np.power(constants.c,2) / np.power(wavelength,5) * 1 / (np.exp(constants.h * constants.c / (wavelength * constants.Boltzmann * T)) - 1)

for k in range(len(T)):
	for j in range(len(wavelengths)):
		B[k,j] = B_fun(wavelengths[j], T[k])

# Plotting
plt.figure(figsize=(7.5,7.5))
plt.title('Blackbody radiation')
plotColours = ['r','g','b','k']
for k in range(len(T)):
	plt.plot(wavelengths*1e6, B[k,:], plotColours[k], label=r'$T = %d$ K' % T[k])

bands = {
	'V' : [0.55, 'green'],
	'I' : [0.8, 'gold'],
	'J' : [1.25, 'darkorange'],
	'H' : [1.635, 'orangered'],
	'K' : [2.2, 'darkred']
}

for key in bands:
	plt.plot([bands[key][0], bands[key][0]], [min(B[0,:]), 1e5], linestyle='dashed', color=bands[key][1])
	plt.text(bands[key][0] + 0.05, 1e-2, key)
# plt.plot([0.55,0.55], [min(B[0,:]), 1e5], linestyle='dashed', color='green', label=r'$V$-band')	
# plt.plot([0.8,0.8], [min(B[0,:]), 1e5], linestyle='dashed', color='gold', label=r'$I$-band')
# plt.plot([1.25,1.25], [min(B[0,:]), 1e5], linestyle='dashed', color='darkorange', label=r'$J$-band')	
# plt.plot([1.635,1.635], [min(B[0,:]), 1e5], linestyle='dashed', color='orangered', label=r'$H$-band')	
# plt.plot([2.2,2.2], [min(B[0,:]), 1e5], linestyle='dashed', color='darkred', label=r'$K$-band')	



plt.xlabel('Wavelength (microns)')
plt.ylabel(r'Spectral radiance (W sr$^{-1}$ m$^{-2}$ m$^{-1}$)')
plt.legend(loc='lower right')
plt.yscale('log')
plt.axis('tight')
plt.show()