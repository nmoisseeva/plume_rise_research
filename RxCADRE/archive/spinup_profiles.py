#testing quality of spinup

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf


#====================INPUT===================
wrfpath = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs_spinup'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
tsample = 30 			#min between samples
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfpath)
wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')

ncdict = wrf.extract_vars(wrfdata, None, ('T','U','V','W','PHB','PH'))
# u = wrf.destagger(ncdict['U'],3)
# v = wrf.destagger(ncdict['V'],2)
# w = wrf.destagger(ncdict['W'],1)

#get geopotential array and convert to height
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
z_ave = np.mean(z, (0,2,3))
nT,nZ,nY,nX = np.shape(z)


plt.figure()
for time in range(0,nT,tsample):
	print time
	plt.plot(ncdict['T'][time,:,int(nY/2),int(nX/2)],z_ave, label='%s min' %time)
plt.plot(ncdict['T'][-1,:,int(nY/2),int(nX/2)],z_ave, label='%s min' %nT)
plt.legend()
plt.xlabel('perturbation of PT [K]')
plt.ylabel('height [m]')
plt.title('Spinup BL Growth')
plt.savefig(fig_dir + 'SpinupBLGrowth.pdf')
plt.show()

