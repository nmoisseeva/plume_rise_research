#script to convert micromet data to sensible heat flux and look at whats happenning
#nmoisseeva@eoas.ubc.ca
#Jan 2019

import numpy as np
import matplotlib.pyplot as plt

datapath = './H.csv'
data = np.genfromtxt(datapath,usecols=(1,2,3,4),skip_header=4,delimiter=',',dtype=float)

freq = 1 #hz
ave_int = 5 #min
num_pts = 60 * freq * ave_int
data_samples = np.shape(data)[0] / num_pts


#samples covering ignition (HARDCODED starting at 12:27pm, with 10am spinup)
bad = range(28, 1000)

H = []

for nSample in range(data_samples):
	if nSample not in bad:
		subset = data[nSample*num_pts:(nSample+1)*num_pts, 2:].T
		min_cov = np.cov(subset)[0,1]
		H.append(min_cov)
	else:
		print 'excluding value'

plt.title('SURFACE HEAT FLUX')
lbl = ['10:00','10:20','10:40','11:00','11:20','11:40','12:00','12:20']
plt.plot(H)
plt.ylim([0,0.3])
plt.xlabel('time (CST)')
plt.ylabel('kinematic heat flux ($\overline{T\'w\'} $) [$K m s^{-1}$]')
ax = plt.gca()
ax.set_xticks(np.arange(0,29,4))
ax.set_xticklabels(lbl)
plt.savefig(fig_dir + 'KinH.pdf')
plt.show()

mean_flux = np.mean(H)
print('Mean flux for pre-ignition period: %s ' %mean_flux)
