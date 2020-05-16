#nmoisseeva@eoas.ubc.class
#sciprt to test spinpup spectra

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import wrf
from scipy.signal import welch
import imp


#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time

#====================INPUT===================


testLvl = 30 			#height level to run the analysis on
#=================end of input===============


print('Extracting NetCDF data from %s ' %plume.wrfdir)
wrfdata = netcdf.netcdf_file(plume.wrfdir + 'wrfout_W5F7R5Tspinup', mode ='r')

ncdict = wrf.extract_vars(wrfdata, None, ('U','V','W','T','PHB','PH'),meta=False)
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)
w = wrf.destagger(ncdict['W'],1)
t = ncdict['T']


print('.....destaggering height data')
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
zstag_ave = np.mean(zstag, (2,3))
z = wrf.destagger(zstag_ave,1)

z_ave = np.mean(z, (0))

nT,nZ1,nY,nX = np.shape(zstag)
nZ = nZ1 - 1


#compare spectra every 2 min
tloop = np.arange(0,nT,2)
u0,v0,w0,t0 = u[0,testLvl,0,0], v[0,testLvl,0,0], w[0,testLvl,0,0],t[0,testLvl,0,0]

plt.figure(figsize=(12,6))
plt.subplot(2,2,1)
for time in tloop:
	uprime = u[time,testLvl,:,:].ravel() - u0
	f,P = welch(uprime, fs=1./plume.dx,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
ulabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(ulabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$u\'$ spetrum')
plt.subplot(2,2,2)
for time in tloop:
	vprime = v[time,testLvl,:,:].ravel() - v0
	f,P = welch(vprime, fs=1./plume.dx,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
vlabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(vlabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$v\'$ spetrum')
plt.subplot(2,2,3)
for time in tloop:
	wprime = w[time,testLvl,:,:].ravel() - w0
	f,P = welch(wprime, fs=1./plume.dx,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
wlabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(wlabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$w\'$ spetrum')
plt.subplot(2,2,4)
for time in tloop:
	tprime = t[time,testLvl,:,:].ravel() - t0
	f,P = welch(tprime, fs=1./plume.dx,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
tlabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(tlabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$t\'$ spetrum')
plt.tight_layout()
# plt.savefig(rx.fig_dir + 'spinup_spectra.pdf')
plt.show()
