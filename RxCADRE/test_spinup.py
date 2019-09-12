# referencing simulation with flight data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf
from scipy.signal import welch
import imp


#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx)		        #force load each time

#====================INPUT===================


lvl = 5 			#height level to run the analysis on
xstep = 40 			#grud spacing in m
#=================end of input===============


print('Extracting NetCDF data from %s ' %rx.spinup_path)
wrfdata = netcdf.netcdf_file(rx.spinup_path, mode ='r')

ncdict = wrf.extract_vars(wrfdata, None, ('U','V','W'))
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)
w = wrf.destagger(ncdict['W'],1)


#get height and destagger
if os.path.isfile(rx.z_path):
    print('.....loading destaggered height array from %s' %rx.z_path)
    z = np.load(rx.z_path)
else:
    print('.....destaggering height data')
    zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
    z = wrf.destagger(zstag,1)
    np.save(rx.z_path, z)

z_ave = np.mean(z, (0,2,3))

nT,nZ,nY,nX = np.shape(z)


#compare spectra every minute
tloop = np.arange(0,nT,120/rx.hist_int)
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
plt.savefig(rx.fig_dir + 'spinup_spectra.pdf')
plt.show()

#compare pre-spinup and post-spinup wind profiles
plt.figure()
plt.plot(u[0,:,0,0], z[0,:,0,0],'r',label='raw profile')
plt.plot(np.mean(np.mean(u[-1,:,:,:],-1),-1), np.mean(np.mean(z[-1,:,:,:],-1),-1),'k', label="post-spinup average profile")
plt.legend()
plt.show()
