# nmoisseeva@eoas.ubc.ca
# January 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf
import cmocean

#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/'
tag = ['W4S400F3R0','W4S400F3R0']
# tag = ['W4S400F3R0']
fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'
lvl = np.arange(0,2500,40)	 			#vertical levels in m
half_width = 20 						#number of grids to use for averaging ALONG fireline
window = [10,65] 						#averaging moving window - number of grids before and after peak flux

#=================end of input===============

print('INTERPOLATION AND AVERAGING SCRIPT FOR PLUME RUNS')
print('===================================')


for nCase,Case in enumerate(tag):
	print('Examining case: %s ' %tag)

	#----------check for interpolated data----------------------------
	interppath = wrfdir + 'interp/wrfinterp_' + tag + '.npy'
	wrfpath = wrfdir + 'wrfout_'+ tag
	wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')

	if os.path.isfile(interppath):
		interpdict = np.load(interppath).item()   # load here the above pickle
		print('Interpolated data found at: %s' %interppath)
		ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX'))
	else:
		print('WARNING: no interpolated data found - generating: SLOW ROUTINE!')

		ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX','W','QVAPOR','T','PHB','PH','U','P','PB','V','XTIME'))

		#get height and destagger vars
		zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
		z = wrf.destagger(zstag,1)
		u = wrf.destagger(ncdict['U'],3)
		v = wrf.destagger(ncdict['V'],2)
		w = wrf.destagger(ncdict['W'],1)
		p = ncdict['P'] + ncdict['PB']

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
					fq = interpolate.interp1d(z_t,ncdict['QVAPOR'][t,:,y,x],fill_value="extrapolate")
					ft = interpolate.interp1d(z_t,ncdict['T'][t,:,y,x],fill_value="extrapolate")
					fw = interpolate.interp1d(z_t,w[t,:,y,x],fill_value="extrapolate")
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

	#convert and average data-------------------------------------------
	ghfx = ncdict['GRNHFX']/1000. 			#convert to kW
	qvapor = interpdict['QVAPOR']*1000.		#convert to g/kg
	temp = interpdict['T']+300. 			#add perturbation and base temperature
	w = interpdict['W']

	#get dimensions
	dimt, dimy, dimx = np.shape(ghfx)
	xsx = int(round(dimy/2.))

	#create fire cross-section averages
	ghfx_mave = np.mean(ghfx[:,xsx-half_width:xsx+half_width,:],1)
	w_mave = np.mean(w[:,:,xsx-half_width:xsx+half_width,:],2)
	#!!!!!!should i look at total or aveage?
	qvapor_mtot = np.nansum(qvapor[:,:,xsx-half_width:xsx+half_width,:],2)

	#create time-average around peak flux--------------------------
	xmax = np.argmax(ghfx_mave,axis=1)
	ghfx_t, qvapor_t, temp_t, w_t, v_t, u_t = [],[],[],[],[],[]
	for nP, pt in enumerate(xmax[1:]):
		ghfx_t.append(ghfx_mave[nP+1,pt-window[0]:pt+window[1]])
		qvapor_t.append(qvapor_mtot[nP+1,:,pt-window[0]:pt+window[1]])
		w_t.append(w_mave[nP+1,:,pt-window[0]:pt+window[1]])
		v_t.append(v_mave[nP+1,:,pt-window[0]:pt+window[1]])
		u_t.append(u_mave[nP+1,:,pt-window[0]:pt+window[1]])
		temp_t.append(temp_mave[nP+1,:,pt-window[0]:pt+window[1]])
	ghfx_tave = np.mean(ghfx_t,0)
	qvapor_tave = np.mean(qvapor_t,0)
	w_tave = np.mean(w_t,0)
	v_tave = np.mean(v_t,0)
	u_tave = np.mean(u_t,0)
	temp_tave = np.mean(temp_t, 0)

	avepath = wrfdir + 'interp/wrfave_' + tag + '.npy'
	avedict = {'GHFX': ghfx_tave, 'T':temp_tave, 'QVAPOR':qvapor_tave, 'W': w_tave, 'U': u_tave, 'V': v_tave}
	np.save(avedict, avepath)
