# Required packages
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

" as a function of r_0 "
# D = np.array([1.0, 2.5, 8.])
# r_0 = np.linspace(1e-2, 20e-2, 100e-2)
# P = np.zeros((len(D), len(r_0)))
# for k in range(len(D)):
# 	P[k,:] = 5.6 * np.exp(-0.1557 * np.power(D[k]/r_0,2))

" Raw numbers "
v_wind = 5		# wind speed at SSO (Smith 2009)
wavelengths = np.array([0.55, 0.8, 1.25, 1.635, 2.2]) * 1e-6
wavelength_I = wavelengths[1]
r_0_I = 12e-2	# best case V-band r_0 (Smith 2009)
tau_0_I = r_0_I / v_wind
D = 1.8			# telescope diameter
P_I = 5.6 * np.exp(-0.1557 * np.power(D / r_0_I, 2))
print "r_0\ttau_0\tP/P_I"
r_0_bands = np.zeros((len(wavelengths)))
for k in range(len(wavelengths)):
	r_0_bands[k] = np.power(wavelengths[k] / wavelength_I, 6./5.) * r_0_I
	tau_0 = r_0_bands[k] / v_wind
	P = 5.6 * np.exp(-0.1557 * np.power(D / r_0_bands[k], 2))
	print "%.2f\t%.2f\t%.10e" % (r_0_bands[k]*1e2, tau_0*1e3, P/P_I)

D = 1.8
r_0 = np.linspace(10e-2, 50e-2, 100)
P = 5.6 * np.exp(-0.1557 * np.power(D/r_0,2))

# Plotting
plt.figure(figsize=(7.5,7.5))
plt.title(r'Lucky image probability as a function of $r_0$')
# for k in range(len(D)):
	# plt.plot(r_0*1e2, P[k,:], label=r'$D = %.1f$ m' % D[k])
plt.plot(r_0*1e2, P, label=r'$D = %.1f$ m' % D)
plt.plot([r_0_I*1e2,r_0_I*1e2], [min(P),1], label=r'$I$-band ($\lambda = 0.800$ $\mu$m)', color='green')
plt.plot([r_0_bands[2]*1e2,r_0_bands[2]*1e2], [min(P),1], label=r'$J$-band ($\lambda = 1.250$ $\mu$m)', color='darkorange')
plt.plot([r_0_bands[3]*1e2,r_0_bands[3]*1e2], [min(P),1], label=r'$H$-band ($\lambda = 1.635$ $\mu$m)', color='orangered')
plt.plot([r_0_bands[4]*1e2,r_0_bands[4]*1e2], [min(P),1], label=r'$K$-band ($\lambda = 2.200$ $\mu$m)', color='darkred')
plt.xlabel(r'$r_0$ (cm)')
plt.ylabel(r'Lucky imaging probability $P$')
plt.ylim([min(P),1])
plt.yscale('log')
plt.xlim([min(r_0)*1e2, max(r_0)*1e2])
# plt.axis('tight')
plt.legend(loc='lower right')
plt.show()

# " as a function of D "
# r_0 = np.array([5.e-2, 10.e-2, 20.e-2, 50.e-2])
# r_0_cm = r_0 * 1e2
# D = np.linspace(0.5, 8, 100)
# P = np.zeros((len(r_0), len(D)))
# for k in range(len(r_0)):
# 	P[k,:] = 5.6 * np.exp(-0.1557 * np.power(D/r_0[k],2))

# # Plotting
# plt.figure(figsize=(7.5,7.5))
# plt.title(r'Lucky image probability as a function of $D$')
# for k in range(len(r_0)):
# 	plt.plot(D, P[k,:], label=r'$r_0 = %d$ cm' % r_0_cm[k])
# plt.xlabel(r'$D$ (m)')
# plt.ylabel(r'Lucky imaging probability $P$')
# plt.ylim([1e-50,1])
# plt.xlim([min(D), max(D)])
# plt.yscale('log')
# plt.legend(loc='lower right')
# plt.show()

# " as a function of D/r_0 "
# ratio = np.linspace(3.5,10,100)
# P = 5.6 * np.exp(-0.1557 * np.power(ratio,2))
# plt.figure(figsize=(7.5,7.5))
# plt.plot(ratio, P, 'r')
# plt.xlabel(r'Ratio $D/r_0$')
# plt.ylabel(r'Lucky imaging probability $P$')
# plt.yscale('log')
# plt.title(r'$P$ vs. $D/r_0$')
# plt.xlim([min(ratio), max(ratio)])
# plt.show()


