# referencing simulation with flight data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf
from scipy.signal import welch


#====================INPUT===================
wrfpath = '/Users/nmoisseeva/data/plume/RxCADRE/wrfout_L2G_spinup_Jan15'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
lvl = 5 			#height level to run the analysis on
tstep = 10 			#timestep in sec
xstep = 40 			#grud spacing in m
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfpath)
wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')

ncdict = wrf.extract_vars(wrfdata, None, ('U','V','W','PHB','PH'))
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)
w = wrf.destagger(ncdict['W'],1)

#get geopotential array and convert to height
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
z_ave = np.mean(z, (0,2,3))
nT,nZ,nY,nX = np.shape(z)

# ui = u[0,:,0,0]
# vi = v[0,:,0,0]
# WS = np.sqrt(ui**2 + vi**2)

#compare spectra every minute
tloop = np.arange(0,nT,120/tstep)
u0,v0,w0 = u[0,lvl,0,0], v[0,lvl,0,0], w[0,lvl,0,0]

plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
for time in tloop:
	uprime = u[time,lvl,:,:].ravel() - u0
	f,P = welch(uprime, fs=1./xstep,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
ulabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(ulabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$u\'$ spetrum')
plt.subplot(1,2,2)
for time in tloop:
	vprime = v[time,lvl,:,:].ravel() - v0
	f,P = welch(vprime, fs=1./xstep,nperseg=320, window='hanning')
	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
plt.gca().set_xticks(f[::30])
vlabels = ['%2.f' %i for i in 1/f[::30]]
plt.gca().set_xticklabels(vlabels)
plt.xlabel('scale [m]')
plt.ylabel('power')
plt.title('$v\'$ spetrum')
plt.tight_layout()
plt.savefig(fig_dir + 'spinup_spectra.pdf')
plt.show()

#compare pre-spinup and post-spinup wind profiles
plt.figure()
plt.plot(u[0,:,0,0], z[0,:,0,0],'r',label='raw profile')
plt.plot(np.mean(np.mean(u[-1,:,:,:],-1),-1), np.mean(np.mean(z[-1,:,:,:],-1),-1),'k', label="post-spinup average profile")
plt.legend()
plt.show()
