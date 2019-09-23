#nmoisseeva@eoas.ubc.ca
#script to compare spinup methods

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import wrf
import imp

#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx)		        #force load each time
#====================INPUT===================
testTime = 140          #min since start
testLvl = 17 			#height level to run the analysis on
blLvl = 28              #level of BL top
#=================end of input===============

#extract data from bubble simulation
print('Extracting NetCDF data from %s ' %rx.spinup_path)
wrfdata = netcdf.netcdf_file(rx.spinup_path, mode ='r')
ncdict = wrf.extract_vars(wrfdata, None, ('W','T','PHB','PH','XTIME'))
time1 = np.where(ncdict['XTIME']==testTime)[0][0]
w1 = ncdict['W'][time1,:blLvl,:,:].ravel()
t1 = (ncdict['T'][time1,:blLvl,:,:] + 300).ravel()
t1BL = (ncdict['T'][time1,:blLvl,:,:] + 300).ravel()

#extract data from perturbation simulation
spinup_tsk = rx.spinup_path[:-6] + 'tsk'
print('Extracting NetCDF data from %s ' %spinup_tsk )
wrfdata_tsk = netcdf.netcdf_file(spinup_tsk, mode ='r')
ncdict_tsk = wrf.extract_vars(wrfdata_tsk, None, ('W','T','PHB','PH','XTIME'))
time2 = np.where(ncdict_tsk['XTIME']==testTime)[0][0]
w2 = ncdict_tsk['W'][time2,:blLvl,:,:].ravel()
t2 = (ncdict_tsk['T'][time2,:blLvl,:,:] + 300).ravel()
t2BL = (ncdict_tsk['T'][time2,:blLvl,:,:] + 300).ravel()

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
t_ave = np.mean(ncdict_tsk['T'][time2,:,:,:],(1,2))
nT,nZ,nY,nX = np.shape(ncdict_tsk['W'])

print('Mean potential temperature (Bubble, TSK): %.3f, %.3f' %(np.mean(t1BL), np.mean(t2BL)))
print('Meadian potential temperature (Bubble, TSK): %.3f, %.3f' %(np.median(t1BL), np.median(t2BL)))
print('25th percentile potential temperature (Bubble, TSK): %.3f, %.3f' %(np.percentile(t1, [25]), np.percentile(t2, [25])))
print('75th percentile potential temperature (Bubble, TSK): %.3f, %.3f' %(np.percentile(t1, [75]), np.percentile(t2, [75])))


plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.title('VERTICAL VELOCITY W')
testW = (w1,w2)
plt.boxplot(testW,showmeans=True,whis=4)
plt.ylabel('w [m/s]')
plt.gca().set_xticklabels(['Bubble','TSK'])
# plt.ylim([-3,6])
plt.subplot(1,2,2)
plt.title('POTENTIAL TEMPERATURE T')
testT = (t1,t2)
plt.boxplot(testT,showmeans=True,whis=4)
plt.ylabel('potential temperature [K]')
plt.gca().set_xticklabels(['Bubble','TSK'])
# plt.ylim([293,294.5])
plt.tight_layout()
plt.savefig(rx.fig_dir + 'SpinupBoxplots.pdf')
plt.show()
#
# #compare spectra every minute
# tloop = np.arange(0,nT,120/rx.hist_int)
# u0,v0,w0 = u[0,testLvl,0,0], v[0,testLvl,0,0], w[0,testLvl,0,0]
#
# plt.figure(figsize=(12,6))
# plt.subplot(1,2,1)
# for time in tloop:
# 	uprime = u[time,testLvl,:,:].ravel() - u0
# 	f,P = welch(uprime, fs=1./xstep,nperseg=320, window='hanning')
# 	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
# plt.gca().set_xticks(f[::30])
# ulabels = ['%2.f' %i for i in 1/f[::30]]
# plt.gca().set_xticklabels(ulabels)
# plt.xlabel('scale [m]')
# plt.ylabel('power')
# plt.title('$u\'$ spetrum')
# plt.subplot(1,2,2)
# for time in tloop:
# 	vprime = v[time,testLvl,:,:].ravel() - v0
# 	f,P = welch(vprime, fs=1./xstep,nperseg=320, window='hanning')
# 	plt.semilogy(f,P, color='k', alpha = float(time)/nT)
# plt.gca().set_xticks(f[::30])
# vlabels = ['%2.f' %i for i in 1/f[::30]]
# plt.gca().set_xticklabels(vlabels)
# plt.xlabel('scale [m]')
# plt.ylabel('power')
# plt.title('$v\'$ spetrum')
# plt.tight_layout()
# plt.savefig(rx.fig_dir + 'spinup_spectra.pdf')
# plt.show()
#
# #compare pre-spinup and post-spinup wind profiles
# plt.figure()
# plt.plot(u[0,:,0,0], z[0,:,0,0],'r',label='raw profile')
# plt.plot(np.mean(np.mean(u[-1,:,:,:],-1),-1), np.mean(np.mean(z[-1,:,:,:],-1),-1),'k', label="post-spinup average profile")
# plt.legend()
# plt.show()
