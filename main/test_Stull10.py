# nmoisseeva@eoad.ubc.ca
# Sept 2017

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import matplotlib.animation as animation
from matplotlib import path 
import os.path
import sys
from matplotlib import animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import wrf
import cmocean

#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/'
tag = 'W4S400F3R0'
# tag = 'W4S400F13R0'
# tag = 'W10S400F3R0'
# tag = 'W4Sn400F3R0'
fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'

dz = 40.
lvl = np.arange(0,2500,dz)	 			#vertical levels in m

#define constants and semi-constants
rho = 1.15 		#average density of smoke plume [kg/m3] (assumes average plume temperature of 30-35C)
Cp = 1000 		#specific heat of of homid smoky air at constant pressure [J/kg]
g = 9.81		#gravity

#=================end of input===============

print('FIRE CROSS-SECTION HEAT FLUX AND W')
print('===================================')


#import data
wrfpath = wrfdir + 'wrfout_'+ tag

print('Extracting NetCDF data from %s ' %wrfpath)
wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')  

#get interpolated WRF data----------------------------------------------------
ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX'))
dx = wrfdata.DX


interppath = wrfdir + 'interp/wrfinterp_' + tag + '.npy'
if os.path.isfile(interppath):
	interpdict = np.load(interppath).item()   # load here the above pickle 
	# qinterp, winterp, interpt = interpdict[()]['QVAPOR'],interpdict[()]['W'], interpdict[()]['T']
	print('Interpolated data found at: %s' %interppath)
else:
	print('ERROR: no interpolated data found - run HFX_W.py to generate!')
	sys.exit()

#convert and average data-------------------------------------------
ghfx = ncdict['GRNHFX']/1000. 				#convert to kW		
qvapor = interpdict['QVAPOR']*1000.		#convert to g/kg		
temp = interpdict['T']+300. 			#add perturbation and base temperature
w = interpdict['W']
u = interpdict['U']

dimt, dimy, dimx = np.shape(ghfx)
xsx = int(round(dimy/2.))
									
#create fire cross-section averages
ghfx_mave = np.mean(ghfx[:,xsx-20:xsx+20,:],1)
w_mave = np.mean(w[:,:,xsx-20:xsx+20,:],2)
u_mave = np.mean(u[:,:,xsx-20:xsx+20,:],2)
qvapor_mtot = np.nansum(qvapor[:,:,xsx-20:xsx+20,:],2)
temp_mave = np.mean(temp[:,:,xsx-20:xsx+20,:],2)


#create time-average around peak flux--------------------------
xmax = np.argmax(ghfx_mave,axis=1)
ghfx_t, qvapor_t, temp_t, w_t, u_t = [],[],[],[],[]
for nP, pt in enumerate(xmax[1:]):
	ghfx_t.append(ghfx_mave[nP+1,pt-10:pt+65])
	qvapor_t.append(qvapor_mtot[nP+1,:,pt-10:pt+65])
	w_t.append(w_mave[nP+1,:,pt-10:pt+65])
	u_t.append(u_mave[nP+1,:,pt-10:pt+65])
	temp_t.append(temp_mave[nP+1,:,pt-10:pt+65])
ghfx_tave = np.mean(ghfx_t,0) 
qvapor_tave = np.mean(qvapor_t,0)
w_tave = np.mean(w_t,0)
u_tave = np.mean(u_t,0)
temp_tave = np.mean(temp_t, 0)

# calculate fireline intensity
fireline = ghfx_tave[ghfx_tave>5] 			#get active fireline (flux > 5000W/m2)
fire_width = len(fireline) * dx 			#get width of the fireline
Iave = fire_width * np.mean(fireline)		#calculate average fireline intensity


#get environmental average U and T - consider bottom layer first
usound = u[0,:,:,:].mean(axis=(1,2))
tsound = temp[0,:,:,:].mean(axis=(1,2))
alpha = fire_width/dz

#test Stull eqn 2
Uf0 = ((2. * g * (alpha**2.) * Iave)/(rho* Cp* tsound[0]))**(1./3.)
u_diff = u_tave[0,:] - usound[0]


#test Stull eqn 10
rhs = (2. * g * alpha**2. ) * Iave/ (rho * Cp * tsound[0]) 
coeffs = [1,usound[0],0,-rhs]
np.roots(coeffs)



#==================================PLOTTING====================================
fig = plt.figure(figsize=(6,6))
ax = plt.gca()
#u contours
im = ax.imshow(u_tave, origin='lower', cmap=cmocean.cm.tempo)
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
plt.savefig(fig_dir + 'U_qv/AVE_hfx_w_qv_%s.pdf' %tag)
plt.show()
print('.....-->saved in: %s' %(fig_dir + 'U_qv/AVE_hfx_w_qv_%s.pdf' %tag))




# #plot contours
# fig = plt.figure(figsize=(6,6))
# ax = plt.gca()
# # ---w contours and colorbar
# cntrf = ax.contourf(w_tave, cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('vertical velocity $[m s^{-2}]$')
# ax.set_xlabel('x grid (#)')
# ax.set_ylabel('height AGL [m]')
# ax.set_yticks(np.arange(0,len(lvl),10))
# ax.set_yticklabels(lvl[::10])
# ax.set_xlim([0,74])
# # ---non-filled vapor contours and colorbar
# cntr = ax.contour(qvapor_tave, cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=1)
# ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
# cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
# cbar.set_label('$H_2O$ mixing ratio $[g kg^{-1}]$',size=8)
# cbar.ax.tick_params(labelsize=8) 
# # ---heat flux
# axh = ax.twinx()
# axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# axh.set_ylim([0,140])
# axh.set_xlim([0,74])
# axh.tick_params(axis='y', colors='red')
# ln = axh.plot(ghfx_tave, 'r-')
# plt.savefig(fig_dir + 'AVE_hfx_w_qv_%s.pdf' %tag)
# plt.show()
# print('.....-->saved in: %s' %(fig_dir + 'AVE_hfx_w_qv_%s.pdf' %tag))


# #plot contours
# fig = plt.figure(figsize=(6,6))
# ax = plt.gca()
# # ---w contours and colorbar
# cntrf = ax.contourf(w_tave, cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('vertical velocity $[m s^{-2}]$')
# ax.set_xlabel('x grid (#)')
# ax.set_ylabel('height AGL [m]')
# ax.set_yticks(np.arange(0,len(lvl),10))
# ax.set_yticklabels(lvl[::10])
# ax.set_xlim([0,74])
# # ---temperature contours and colorbar
# cntr = ax.contour(temp_tave, cmap=plt.cm.YlOrRd,levels=np.arange(300,306.1,0.5),linewidths=1)
# ains = inset_axes(plt.gca(), width='50%', height='2%', loc=1)
# cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
# cbar.set_label('potential temperature $[K]$',size=8)
# cbar.ax.tick_params(labelsize=8) 
# # ---heat flux
# axh = ax.twinx()
# axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# axh.set_ylim([0,140])
# axh.set_xlim([0,74])
# axh.tick_params(axis='y', colors='red')
# ln = axh.plot(ghfx_tave, 'r-')
# plt.savefig(fig_dir + 'AVE_hfx_w_temp_%s.pdf' %tag)
# plt.show()
# print('.....-->saved in: %s' %(fig_dir + 'AVE_hfx_w_temp_%s.pdf' %tag))



# #create animation of the HFX and W ------------------------------------
# print('.....creating crossection of HFX and W animation')
# fig = plt.figure()
# plt.title('AVE CROSS-SECTION HFX and W')
# ax1 = plt.gca()
# # create initial frame
# ln1 = ax1.plot(ghfx_mave[0,:], 'r-')
# ax1.set_xlabel('x grid (#)')
# ax1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# ax1.set_ylim([0,200])
# ax2 = ax1.twinx()
# ln2 = ax2.plot(w_mave[0,5,:], 'b:')
# ax2.set_ylabel('vertical velocity $[m s^{-2}]$', color='b')
# ax2.set_ylim([0,3])

# def update_plot(n, ghfx,w):
# 	ax1.clear()
# 	ax1.set_xlabel('x grid (#)')
# 	ax1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# 	ax1.set_ylim([0,200])
# 	ax2.clear()
# 	ax2.set_ylabel('vertical velocity $[m s^{-2}]$', color='b')
# 	ax2.set_ylim([0,5])
# 	ln1 = ax1.plot(ghfx_mave[n,:], 'r-')
# 	ln2 = ax2.plot(w_mave[n,3,:], 'b:')
# 	# time_text.set_text('Time (sec) = %s' %(n*dt))
# 	return ln1, ln2,

# #plot all frames
# ani=animation.FuncAnimation(fig, update_plot, dimt, fargs=(ghfx_mave,w_mave), interval=3)
# ani.save(fig_dir + 'CW_hfx_w_%s.mp4' %tag, writer='ffmpeg',fps=10, dpi=200)
# print('.....-->saved in: %s' %(fig_dir + 'CW_hfx_w_%s.mp4' %tag))



# #create animation of vertical velocity contours and water vapor-----------------------
# print('.....creating vertical crossection of W + H2O animation')
# fig = plt.figure(figsize=(6,6))
# ax = plt.gca()
# # create initial frame
# # ---w contours and colorbar
# cntrf = ax.contourf(w_mave[0,:,:], cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5))
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('vertical velocity $[m s^{-2}]$')
# ax.set_xlabel('x grid (#)')
# ax.set_ylabel('height AGL [m]')
# ax.set_xlim([0,dimx])
# # ---non-filled vapor contours and colorbar
# cntr = ax.contour(qvapor_mtot[0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
# ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
# cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
# cbar.set_label('$H_2O$ mixing ratio $[g kg^{-1}]$',size=8)
# cbar.ax.tick_params(labelsize=8) 
# # ---heat flux
# axh = ax.twinx()
# axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# axh.set_ylim([0,140])
# axh.set_xlim([0,dimx])
# axh.tick_params(axis='y', colors='red')
# ln = axh.plot(ghfx_mave[0,:], 'r-')

# def update_plot(n,w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr):
#         ax.clear()
#         cntrf = ax.contourf(w_mave[n,:,:],cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
#         cntr = ax.contour(qvapor_mtot[n,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=0.6)
#         ax.set_xlabel('x grid (#)')
#         ax.set_ylabel('height AGL [m]')
#         ax.set_yticks(np.arange(0,len(lvl),10))
#         ax.set_yticklabels(lvl[::10])
#         axh.clear()
#         axh.set_ylim([0,140])
#         axh.set_xlim([0,dimx])
#         ln = axh.plot(ghfx_mave[n,:], 'r-')
#         axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
#         # time_text.set_text('Time (sec) = %s' %(n*dt))
#         return cntrf, ln, cntr,

# #plot all frames
# ani=animation.FuncAnimation(fig, update_plot, dimt, fargs=(w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr), interval=3)
# ani.save(fig_dir + 'Vert_w_qvapor_%s.mp4' %tag, writer='ffmpeg',fps=10, dpi=250)
# plt.close()
# print('.....saved in: %s' %(fig_dir + 'Vert_w_qvapor_%s.mp4' %tag))


# #create animation of vertical velocity contours and temperature--------------------
# print('.....creating vertical crossection of W + temperature animation')
# fig = plt.figure(figsize=(6,6))
# ax = plt.gca()
# #create initial frame
# # ---w contours and colorbar
# cntrf = ax.contourf(w_mave[0,:,:], cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('vertical velocity $[m s^{-2}]$')
# ax.set_xlabel('x grid no.')
# ax.set_ylabel('height AGL [m]')
# # ---temperature contours and colorbar
# cntr = ax.contour(temp_mave[0,:,:], cmap=plt.cm.YlOrRd,levels=np.arange(300,305.1,0.5),linewidths=4)
# ains = inset_axes(plt.gca(), width='50%', height='2%', loc=1)
# cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
# cbar.set_label('potential temperature $[K]$',size=8)
# cbar.ax.tick_params(labelsize=8) 
# #---heat flux
# axh = ax.twinx()
# axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# axh.set_ylim([0,140])
# axh.set_xlim([0,dimx])
# axh.tick_params(axis='y', colors='red')
# ln = axh.plot(ghfx_mave[0,:], 'r-')

# def update_plot(n,w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr):
# 	ax.clear()
# 	cntrf = ax.contourf(w_mave[n,:,:], cmap=plt.cm.PRGn_r, levels=np.arange(-7,7.1,0.5),extend='both')
# 	cntr = ax.contour(temp_mave[n,:,:], cmap=plt.cm.YlOrRd,levels=np.arange(300,305.1,0.5),linewidths=0.5) 	#plume temperature
# 	ax.set_xlabel('x grid (#)')
# 	ax.set_ylabel('height AGL [m]')
# 	ax.set_yticks(np.arange(0,len(lvl),10))
# 	ax.set_yticklabels(lvl[::10])
# 	axh.clear()
# 	axh.set_ylim([0,140])
# 	axh.set_xlim([0,dimx])
# 	ln = axh.plot(ghfx_mave[n,:], 'r-')
# 	axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# 	# time_text.set_text('Time (sec) = %s' %(n*dt))
# 	return cntrf, ln, cntr,


# #plot all frames
# ani=animation.FuncAnimation(fig, update_plot, dimt, fargs=(w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr), interval=3)
# ani.save(fig_dir + 'Vert_w_temp_%s.mp4' %tag, writer='ffmpeg',fps=10, dpi=250)
# plt.close()
# print('.....saved in: %s' %(fig_dir + 'Vert_w_temp_%s.mp4' %tag))



# #create animation of plume contours and temperature--------------------
# print('.....creating vertical crossection of temperature + plume animation')
# fig = plt.figure(figsize=(6,6))
# ax = plt.gca()
# # create initial frame
# # ---w contours and colorbar
# cntrf = ax.contourf(temp_mave[0,:,:], cmap=cmocean.cm.thermal, levels=np.arange(300,316,1.),extend='both')
# cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
# cbarf.set_label('potential temperature $[K]$')
# ax.set_xlabel('x grid no.')
# ax.set_ylabel('height AGL [m]')
# # ---temperature contours and colorbar
# cntr = ax.contour(qvapor_mtot[0,:,:], cmap=cmocean.cm.gray_r,levels=np.arange(0,2.1,0.3),linewidths=4)
# ains = inset_axes(plt.gca(), width='50%', height='2%', loc=1)
# cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
# cbar.set_label('vertical velocity $[m s^{-2}]$',size=8)
# cbar.ax.tick_params(labelsize=8) 
# # ---heat flux
# axh = ax.twinx()
# axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# axh.set_ylim([0,140])
# axh.set_xlim([0,dimx])
# axh.tick_params(axis='y', colors='red')
# ln = axh.plot(ghfx_mave[0,:], 'r-')

# def update_plot(n,w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr):
# 	ax.clear()
# 	cntrf = ax.contourf(temp_mave[n,:,:],cmap=cmocean.cm.thermal,levels=np.arange(300,316,1.),extend='both')
# 	cntr = ax.contour(qvapor_mtot[n,:,:], cmap=cmocean.cm.gray_r,levels=np.arange(0,2.1,0.3),linewidths=0.5) 	#plume temperature
# 	ax.set_xlabel('x grid (#)')
# 	ax.set_ylabel('height AGL [m]')
# 	ax.set_yticks(np.arange(0,len(lvl),10))
# 	ax.set_yticklabels(lvl[::10])
# 	axh.clear()
# 	axh.set_ylim([0,140])
# 	axh.set_xlim([0,dimx])
# 	ln = axh.plot(ghfx_mave[n,:], 'r-')
# 	axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
# 	# time_text.set_text('Time (sec) = %s' %(n*dt))
# 	return cntrf, ln, cntr,


# #plot all frames
# ani=animation.FuncAnimation(fig, update_plot, dimt, fargs=(w_mave,cntrf,ghfx_mave,qvapor_mtot,cntr), interval=3)
# ani.save(fig_dir + 'Vert_temp_qv_%s.mp4' %tag, writer='ffmpeg',fps=10, dpi=250)
# plt.close()
# print('.....saved in: %s' %(fig_dir + 'Vert_temp_qv_%s.mp4' %tag))

