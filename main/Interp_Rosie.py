# nmoisseeva@eoas.ubc.ca
# Sept 2017

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
from matplotlib import path
import os.path
import wrf

#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/'
tag = 'W4S400F3R0'
lvl = np.arange(0,2500,40)	 			#vertical levels in m

#=================end of input===============

print('FIRE CROSS-SECTION HEAT FLUX AND W')
print('===================================')

#import data
wrfpath = wrfdir + 'wrfout_'+ tag

print('Extracting NetCDF data from %s ' %wrfpath)
wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')

#prep WRF data----------------------------------------------------
ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX','W','QVAPOR','T','PHB','PH','U','P','PB','V'))

#get height and destagger vars
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)
p = ncdict['P'] + ncdict['PB']

interppath = wrfdir + 'interp/wrfinterp_' + tag + '.npy'
if os.path.isfile(interppath):
	interpdict = np.load(interppath).item()
	print('Interpolated data found at: %s' %interppath)
else:
	print('WARNING: no interpolated data found - generating: SLOW ROUTINE!')
	nT,nZ,nY,nX = np.shape(z)
	qinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	winterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	uinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	vinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	tinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	pinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	for t in range(nT):
		print('.... tsetp = %s/%s' %(t,nT))
		for y in range(nY):
			for x in range(nX):
				z_t = z[t,:,y,x]
				zstag_t = zstag[t,:,y,x]
				fq = interpolate.interp1d(z_t,ncdict['QVAPOR'][t,:,y,x],fill_value="extrapolate")
				fw = interpolate.interp1d(zstag_t,ncdict['W'][t,:,y,x],fill_value="extrapolate")
				ft = interpolate.interp1d(z_t,ncdict['T'][t,:,y,x],fill_value="extrapolate")
				fu = interpolate.interp1d(z_t,u[t,:,y,x],fill_value="extrapolate")
				fv = interpolate.interp1d(z_t,v[t,:,y,x],fill_value="extrapolate")
				fp = interpolate.interp1d(z_t,p[t,:,y,x],fill_value="extrapolate")
				qinterp[t,:,y,x] = fq(lvl)
				winterp[t,:,y,x] = fw(lvl)
				tinterp[t,:,y,x] = ft(lvl)
				uinterp[t,:,y,x] = fu(lvl)
				vinterp[t,:,y,x] = fv(lvl)
				pinterp[t,:,y,x] = fp(lvl)
	interpdict = {'QVAPOR': qinterp, 'W':winterp, 'T':tinterp, 'U':uinterp,'P':pinterp, 'V':vinterp}
	np.save(interppath, interpdict)
	print('Interpolated data saved as: %s' %interppath)
