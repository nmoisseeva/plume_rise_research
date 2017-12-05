# nmoisseeva@eoad.ubc.ca
# Dec 2017

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
from matplotlib import path 
import os.path
import wrf


#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/' 	#directory on local Mac
# run = ['W4S400F3R0', 'W4S400F13R0','W10S400F3R0','W4Sn400F3R0']
run = ['W4S400F3R0', 'W4S400F13R0']

interpdir = 'interp/wrfinterp_'
lvl = np.arange(0,2500,40)	 			#vertical levels in m

#=================end of input===============

print('BULK INTERPOLATION OF LES PLUME DATA')
print('======================================')


for iTag, tag in enumerate(run):
	#import data
	wrfpath = wrfdir + 'wrfout_'+ tag

	print('-------------Run: %s--------------' %tag)
	print('.... Extracting data from raw netcdf')
	wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')  

	#prep WRF data----------------------------------------------------
	ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX','W','QVAPOR','T','PHB','PH','U','V','P','PB'))

	#get height and destagger vars
	zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
	z = wrf.destagger(zstag,1)
	u = wrf.destagger(ncdict['U'],3)
	v = wrf.destagger(ncdict['V'],2)
	p = ncdict['P'] + ncdict['PB']

	interppath = wrfdir + interpdir + tag + '.npy'
	
	print('.... Destaggering and interpolating data')
	nT,nZ,nY,nX = np.shape(z)
	qinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	winterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	uinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	vinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	tinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	pinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	for t in range(nT):
		print('....... tsetp = %s/%s' %(t,nT))
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
	interpdict = {'QVAPOR': qinterp, 'W':winterp, 'T':tinterp, 'U':uinterp,'V':vinterp,'P':pinterp}
	np.save(interppath, interpdict)
	print('Interpolated data saved as: %s' %interppath)


