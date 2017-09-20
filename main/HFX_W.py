# nmoisseeva@eoad.ubc.ca
# Sept 2017

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
# from scipy import interpolate
# from scipy.ndimage.interpolation import rotate
# from scipy.spatial import KDTree
import matplotlib.animation as animation
from matplotlib import path 
# from mpl_toolkits import basemap
# import mpl_toolkits.basemap.pyproj as pyproj
import os.path
# import pickle
# import mpl_toolkits.mplot3d.axes3d as p3
# import mpl_toolkits.mplot3d as a3
from matplotlib import animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#====================INPUT===================
wrfdata = '/Users/nmoisseeva/data/plume/main/wrfout_W4S400F3R0'
fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'


#=================end of input===============

print('FIRE CROSS-SECTION HEAT FLUX AND W')
print('===================================')


#import data
print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
ghfx = np.copy(nc_data.variables['GRNHFX'][:,:,:])		#ground heat flux
w = np.copy(nc_data.variables['W'][:,:,:,:]) 			#vertical velocity 
qvapor = np.copy(nc_data.variables['QVAPOR'][:,:,:,:])	#water vapor
temp = np.copy(nc_data.variables['T'][:,:,:,:])			#temperature
xtime = xtime = nc_data.variables['XTIME'][:] 			#time in minutes

#get cross-section index and concert to kW
dimt, dimy, dimx = np.shape(ghfx)
xsx = int(round(dimy/2.))
ghfx = ghfx/1000.										#convert to kW
qvapor = qvapor*1000. 									#convert to g/kg
temp = temp + 300. 										#convert to total temperature

#creat fire cross-section averages
ghfx_mave = np.mean(ghfx[:,xsx-20:xsx+20,:],1)
w_mave = np.mean(w[:,:,xsx-20:xsx+20,:],2)
qvapor_mtot = np.nansum(qvapor[:,:,xsx-20:xsx+20,:],2)
temp_mave = np.mean(temp[:,:,xsx-20:xsx+20,:],2)


# #-----------------VERTICAL INTERPOLATION-----------------!!!!!NOT DONE
# numLvl = len(lvl)

# #open/generate vertically interpolated data
# if os.path.isfile(interp_path):
# 	qinterp = np.load(interp_path)   # load here the above pickle 
# 	print('Interpolated data found at: %s' %interp_path)
# else:
# 	print('WARNING: no interpolated vapour data found - generating: SLOW ROUTINE!')
# 	qvcopy = np.copy(nc_data.variables['QVAPOR'][:,:,:,:])
# 	qinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
# 	for t in range(nT):
# 		print('.... tsetp = %s/%s' %(t,nT))
# 		for y in range(nY):
# 			for x in range(nX):
# 				tempz = z[t,:,y,x]
# 				interpz = (tempz[:-1]+tempz[1:])/2.
# 				f = interpolate.interp1d(interpz,qvcopy[t,:,y,x],fill_value="extrapolate")
# 				qinterp[t,:,y,x] = f(lvl)
# 	np.save(interp_path, qinterp)
# 	print('Interpolated data saved as: %s' %interp_path)


#create animation of the HFX and W ------------------------------------
print('.....creating crossection of HFX and W animation')
fig = plt.figure()
plt.title('AVE CROSS-SECTION HFX and W')
ax1 = plt.gca()
# create initial frame
ln1 = ax1.plot(ghfx_mave[0,:], 'r-')
ax1.set_xlabel('x grid (#)')
ax1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
ax1.set_ylim([0,150])
ax2 = ax1.twinx()
ln2 = ax2.plot(w_mave[0,2,:], 'b:')
ax2.set_ylabel('vertical velocity $[m s^{-2}]$', color='b')
ax2.set_ylim([0,2])

def update_plot(n, ghfx,w):
	ax1.clear()
	ax1.set_xlabel('x grid (#)')
	ax1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
	ax1.set_ylim([0,150])
	ax2.clear()
	ax2.set_ylabel('vertical velocity $[m s^{-2}]$', color='b')
	ax2.set_ylim([0,2])
	ln1 = ax1.plot(ghfx_mave[n,:], 'r-')
	ln2 = ax2.plot(w_mave[n,2,:], 'b:')
	# time_text.set_text('Time (sec) = %s' %(n*dt))
	return ln1, ln2,

#plot all frames
ani=animation.FuncAnimation(fig, update_plot, len(xtime), fargs=(ghfx,w), interval=3)
ani.save(fig_dir + 'CW_hfx_w.mp4', writer='ffmpeg',fps=10, dpi=200)
plt.show()
print('.....-->saved in: %s' %(fig_dir + 'CW_hfx_w.mp4'))




#create animation of vertical velocity contours and water vapor-----------------------
print('.....creating vertical crossection of W + H2O animation')
fig = plt.figure(figsize=(6,6))
ax = plt.gca()
# create initial frame
# ---w contours and colorbar
cntrf = ax.contourf(w_mave[0,:,:], cmap=plt.cm.PRGn_r, levels=np.arange(-5,5.1,0.5))
cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
cbarf.set_label('vertical velocity $[m s^{-2}]$')
ax.set_xlabel('x grid (#)')
ax.set_xlim([0,dimx])
# ---non-filled vapor contours and colorbar
cntr = ax.contour(qvapor_mtot[0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
cbar.set_label('$H_2O$ mixing ratio $[g kg^{-1}]$',size=8)
cbar.ax.tick_params(labelsize=8) 
# ---heat flux
axh = ax.twinx()
axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
axh.set_ylim([0,150])
axh.set_xlim([0,dimx])
axh.tick_params(axis='y', colors='red')
ln = axh.plot(ghfx_mave[0,:], 'r-')
plt.tight_layout()

def update_plot(n,w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr):
	ax.clear()
	cntrf = ax.contourf(w_mave[n,:,:],cmap=plt.cm.PRGn_r, levels=np.arange(-5,5.1,0.5),extend='both')
	cntr = ax.contour(qvapor_mtot[n,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=0.6) 	
	ax.set_xlabel('x grid (#)')
	axh.clear()
	axh.set_ylim([0,150])
	axh.set_xlim([0,dimx])
	ln = axh.plot(ghfx_mave[n,:], 'r-')
	axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
	# time_text.set_text('Time (sec) = %s' %(n*dt))
	return cntrf, ln, cntr,

#plot all frames
ani=animation.FuncAnimation(fig, update_plot, len(xtime), fargs=(w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr), interval=3)
ani.save(fig_dir + 'Vert_w_qvapor.mp4', writer='ffmpeg',fps=10, dpi=250)
plt.close()
print('.....saved in: %s' %(fig_dir + 'Vert_w_qvapor.mp4'))


#create animation of vertical velocity contours and temperature--------------------
print('.....creating vertical crossection of W + temperature animation')
fig = plt.figure(figsize=(6,6))
ax = plt.gca()
# create initial frame
# ---w contours and colorbar
cntrf = ax.contourf(w_mave[0,:,:], cmap=plt.cm.PRGn_r, levels=np.arange(-5,5.1,0.5),extend='both')
cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
cbarf.set_label('vertical velocity $[m s^{-2}]$')
ax.set_xlabel('x grid (#)')
ax.set_xlim([0,dimx])
# ---temperature contours and colorbar
cntr = ax.contour(qvapor_mtot[0,:,:], cmap=plt.cm.YlOrRd,levels=np.arange(300,305.1,0.5),linewidths=4)
ains = inset_axes(plt.gca(), width='50%', height='2%', loc=1)
cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
cbar.set_label('potential temperature $[K]$',size=8)
cbar.ax.tick_params(labelsize=8) 
# ---heat flux
axh = ax.twinx()
axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
axh.set_ylim([0,150])
axh.set_xlim([0,dimx])
axh.tick_params(axis='y', colors='red')
ln = axh.plot(ghfx_mave[0,:], 'r-')
plt.tight_layout()

def update_plot(n,w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr):
	ax.clear()
	cntrf = ax.contourf(w_mave[n,:,:],cmap=plt.cm.PRGn_r, levels=np.arange(-5,5.1,0.5),extend='both')
	cntr = ax.contour(temp_mave[n,:,:], cmap=plt.cm.YlOrRd,levels=np.arange(300,305.1,0.5),linewidths=0.5) 	#plume temperature
	ax.set_xlabel('x grid (#)')
	axh.clear()
	axh.set_ylim([0,150])
	axh.set_xlim([0,dimx])
	ln = axh.plot(ghfx_mave[n,:], 'r-')
	axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
	# time_text.set_text('Time (sec) = %s' %(n*dt))
	return cntrf, ln, cntr,

#plot all frames
ani=animation.FuncAnimation(fig, update_plot, len(xtime), fargs=(w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr), interval=3)
ani.save(fig_dir + 'Vert_w_temp.mp4', writer='ffmpeg',fps=10, dpi=250)
plt.close()
print('.....saved in: %s' %(fig_dir + 'Vert_w_temp.mp4'))
