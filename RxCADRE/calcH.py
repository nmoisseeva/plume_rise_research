#script to convert micromet data to sensible heat flux and look at whats happenning
#nmoisseeva@eoas.ubc.ca
#Jan 2019

import numpy as np
import matplotlib.pyplot as plt

datapath = './U.csv'
data = np.genfromtxt(datapath,usecols=(1,2,3,4,5,6,7,8,9),skip_header=3,delimiter=',',dtype=float)


freq = 30 #hz
ave_int = 1 #min
num_pts = 60 * freq * ave_int
data_samples = np.shape(data)[0] / num_pts


#samples covering ignition
# bad = range(160,190)
# bad = []
bad = range(14, 20)

H = []

for nSample in range(data_samples):
	if nSample not in bad:
		subset = data[nSample*num_pts:(nSample+1)*num_pts, 2:].T
		min_cov = np.cov(subset)[0,1]
		H.append(min_cov)
	else:
		print 'excluding value'

plt.plot(H)
plt.show()

mean_flux = np.mean(H)
print('Mean flux for pre-ignition period: %s ' %mean_flux)
