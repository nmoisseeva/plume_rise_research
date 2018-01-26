# nmoisseeva@eoas.ubc.ca
# January 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
# from scipy.ndimage.interpolation import rotate
# from scipy.spatial import KDTree
# import matplotlib.animation as animation
from matplotlib import path
# from mpl_toolkits import basemap
# import mpl_toolkits.basemap.pyproj as pyproj
import os.path
# import pickle
# import mpl_toolkits.mplot3d.axes3d as p3
# import mpl_toolkits.mplot3d as a3
# from matplotlib import animation
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import wrf
import cmocean

#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/'
tag = 'W4S400F3R0'
# tag = 'W4S400F13R0'
# tag = 'W10S400F3R0'
# tag = 'W4Sn400F3R0'
fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'

lvl = np.arange(0,2500,40)	 			#vertical levels in m

#=================end of input===============

print('ANALYSIS OF PLUME CENTERLINE')
print('===================================')
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

#convert and average data-------------------------------------------
ghfx = ncdict['GRNHFX']/1000. 			#convert to kW
qvapor = interpdict['QVAPOR']*1000.		#convert to g/kg
temp = interpdict['T']+300. 			#add perturbation and base temperature
w = interpdict['W']
# v = interpdict['V']

dimt, dimy, dimx = np.shape(ghfx)
xsx = int(round(dimy/2.))

#create fire cross-section averages
ghfx_mave = np.mean(ghfx[:,xsx-20:xsx+20,:],1)
w_mave = np.mean(w[:,:,xsx-20:xsx+20,:],2)
# temp_mave = np.mean(temp[:,:,xsx-20:xsx+20,:],2)
#!!!!!!should i look at total or aveage?
qvapor_mtot = np.nansum(qvapor[:,:,xsx-20:xsx+20,:],2)
# v_mave = np.mean(v[:,:,xsx:xsx+20:,:],2)
# u_mave = np.mean(w[:,:,xsx-20:xsx+20,:],2)

#create time-average around peak flux--------------------------
xmax = np.argmax(ghfx_mave,axis=1)
ghfx_t, qvapor_t, temp_t, w_t, v_t, u_t = [],[],[],[],[],[]
for nP, pt in enumerate(xmax[1:]):
	ghfx_t.append(ghfx_mave[nP+1,pt-10:pt+65])
	qvapor_t.append(qvapor_mtot[nP+1,:,pt-10:pt+65])
	w_t.append(w_mave[nP+1,:,pt-10:pt+65])
	# v_t.append(v_mave[nP+1,:,pt-10:pt+65])
	# u_t.append(u_mave[nP+1,:,pt-10:pt+65])
	# temp_t.append(temp_mave[nP+1,:,pt-10:pt+65])
ghfx_tave = np.mean(ghfx_t,0)
qvapor_tave = np.mean(qvapor_t,0)
w_tave = np.mean(w_t,0)
# v_tave = np.mean(v_t,0)
# u_tave = np.mean(u_t,0)
# temp_tave = np.mean(temp_t, 0)

#extract profiles of max w and q
wave_plume = w_tave.copy()
wave_plume[qvapor_tave<0.1] = np.nan 		#mask where there is no plume
wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
qmax_profile = np.nanmax(qvapor_tave,1) 	#get max q profile
qmax_idx = np.nanargmax(qvapor_tave[np.isfinite(qmax_profile)],1)		#get donwind location

#create a vertical profile based on downwind location of maximum vertical lift and max q
qslicew_idx = np.max(wmax_idx)
qslicew = qvapor_tave[:,qslicew_idx]
qsliceq_idx = np.max(qmax_idx)
qsliceq = qvapor_tave[:,qsliceq_idx]

#===========================plotting===========================
#vertical concentration slice at donwind locations of wmax and qmax
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.title('DOWNIND LOCATION OF MAX PROFILES')
plt.plot(wmax_idx,lvl[np.isfinite(wmax_profile)],label='$w_{max}$')
plt.plot(qmax_idx,lvl[np.isfinite(qmax_profile)],'--',label='$q_{max}$')
plt.xlabel('x-grid [#]')
plt.ylabel('height [m]')
plt.legend()
plt.subplot(1,2,2)
plt.title('Q CONCENTRATION DOWNWIND (SLICE)')
plt.plot(qslicew, lvl, label='based on $w_{max}$')
plt.plot(qsliceq, lvl, '--', label='based on $q_{max}$')
plt.xlabel('water vapor [g/kg]')
plt.ylabel('height [m]')
plt.legend()
plt.show()
plt.close()


#plot vertical velocity profile 


#plot velocity distribution filled
	#plume contours
	#scatter points of the wmax and qmax cernterlines - check that it's mid-plume, why do they differ?

#==================================PLOTTING====================================
#plot contours
fig = plt.figure(figsize=(6,6))
ax = plt.gca()
# ---w contours and colorbar
# cntrf = ax.contourf(w_tave, cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('vertical velocity $[m s^{-2}]$')
im = ax.imshow(w_tave, origin='lower', cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
cbari.set_label('horizontal velocity u $[m s^{-2}]$')
ax.set_xlabel('x grid (#)')
ax.set_ylabel('height AGL [m]')
ax.set_yticks(np.arange(0,len(lvl),10))
ax.set_yticklabels(lvl[::10])
ax.set_xlim([0,74])
# ---non-filled vapor contours and colorbar
cntr = ax.contour(qvapor_tave, cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=1)
ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
cbar.set_label('$H_2O$ mixing ratio $[g kg^{-1}]$',size=8)
cbar.ax.tick_params(labelsize=8)
# ---heat flux
axh = ax.twinx()
axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
axh.set_ylim([0,140])
axh.set_xlim([0,74])
axh.tick_params(axis='y', colors='red')
ln = axh.plot(ghfx_tave, 'r-')
plt.savefig(fig_dir + 'TimeAverage/AVE_hfx_w_qv_%s.pdf' %tag)
plt.show()
print('.....-->saved in: %s' %(fig_dir + 'TimeAverage/AVE_hfx_w_qv_%s.pdf' %tag))


#plot contours
fig = plt.figure(figsize=(6,6))
ax = plt.gca()
# ---w contours and colorbar
cntrf = ax.contourf(w_tave, cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
cbarf.set_label('vertical velocity $[m s^{-2}]$')
ax.set_xlabel('x grid (#)')
ax.set_ylabel('height AGL [m]')
ax.set_yticks(np.arange(0,len(lvl),10))
ax.set_yticklabels(lvl[::10])
ax.set_xlim([0,74])
# ---temperature contours and colorbar
cntr = ax.contour(temp_tave, cmap=plt.cm.YlOrRd,levels=np.arange(300,306.1,0.5),linewidths=1)
ains = inset_axes(plt.gca(), width='50%', height='2%', loc=1)
cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
cbar.set_label('potential temperature $[K]$',size=8)
cbar.ax.tick_params(labelsize=8)
# ---heat flux
axh = ax.twinx()
axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
axh.set_ylim([0,140])
axh.set_xlim([0,74])
axh.tick_params(axis='y', colors='red')
ln = axh.plot(ghfx_tave, 'r-')
plt.savefig(fig_dir + 'TimeAverage/AVE_hfx_w_temp_%s.pdf' %tag)
plt.show()
print('.....-->saved in: %s' %(fig_dir + 'TimeAverage/AVE_hfx_w_temp_%s.pdf' %tag))
