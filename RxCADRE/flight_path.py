# referencing simulation with flight data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
from scipy.ndimage.interpolation import rotate
from scipy.spatial import KDTree
import matplotlib.animation as animation
from matplotlib import path
from mpl_toolkits import basemap
import mpl_toolkits.basemap.pyproj as pyproj
import os.path
import pickle
import mpl_toolkits.mplot3d.axes3d as p3
import mpl_toolkits.mplot3d as a3
from matplotlib import animation



#====================INPUT===================
wrfdata = '/Users/nmoisseeva/data/plume/RxCADRE/regrid/wrfout_01-02-2018_regrid'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
bounds_shape = '/Users/nmoisseeva/data/qgis/LG2012_WGS'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
emis_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/Emissions_L2G_20121110.csv'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_01-02-2018_interp.npy'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile

# ll_utm = np.array([519500,3377000])		#lower left corner of the domain in utm
ll_utm = np.array([518300,3377000]) 	#Jan 2018
basemap_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])

lvl = np.arange(0,1700,20) 				#
emis_excl = 0 							#number of samples to excluded (from the END!)
sfc_hgt = 62 							#surface height MSL (m)
# runstart = '12:27:00' 				#start time (if restart run time of inital simulation)
runstart = '12:00:00' 					#start time (if restart run time of inital simulation)
runend = '13:00:00'
corskcrew_ssm= [47183,47521]			#start and end of corskcrew maneuver in ssm
# bg_cork_ssm = [43975,44371] 			#start and end time of pre-burn corkscrew for background
garage_ssm = [45237,47133] 				#start and end time of garage profile
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')

#get geopotential array and convrt to height
z = (nc_data.variables['PHB'][:,:,:,:] + nc_data.variables['PH'][:,:,:,:]) / 9.81

#get domain dimensions
WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]
nT,nZ,nY,nX = np.shape(z)

#open/generate basemap
if os.path.isfile(basemap_path):
	bm = pickle.load(open(basemap_path,'rb'))   # load here the above pickle
	print('Domain basemap found at: %s' %basemap_path)
else:
	print('WARNING: no existing basemaps found: configuring a new basemap')
	bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)
	pickle.dump(bm,open(basemap_path,'wb'),-1)  	# pickle the new map for later
	print('.....New basemap instance saved as: %s' %basemap_path)

# Sanity check: import shape file
# polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)
# fireim = nc_data.variables['GRNHFX'][25,:,:]
# bm.imshow(fireim)
# plt.show()

#extract model time info
tsec = nc_data.variables['XTIME'][:] * 60. 		#get time in seconds since run start
model_ssm = int(runstart[0:2])*3600 + int(runstart[3:5])*60

#==========================VERTICAL INTERPOLATION============================
numLvl = len(lvl)

#open/generate vertically interpolated data
if os.path.isfile(interp_path):
	qinterp = np.load(interp_path)   # load here the above pickle
	print('Interpolated data found at: %s' %interp_path)
else:
	print('WARNING: no interpolated vapour data found - generating: SLOW ROUTINE!')
	qvcopy = np.copy(nc_data.variables['QVAPOR'][:,:,:,:])
	qinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
	for t in range(nT):
		print('.... tsetp = %s/%s' %(t,nT))
		for y in range(nY):
			for x in range(nX):
				tempz = z[t,:,y,x]
				interpz = (tempz[:-1]+tempz[1:])/2.
				f = interpolate.interp1d(interpz,qvcopy[t,:,y,x],fill_value="extrapolate")
				qinterp[t,:,y,x] = f(lvl)
	np.save(interp_path, qinterp)
	print('Interpolated data saved as: %s' %interp_path)

#================================DISPERSION==================================
#extract and format dispersion data
print('Importing dispersion data from %s' %disp_data)
disp_dict = {}
disp_array = np.genfromtxt(disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')

start_idx = np.argmin(abs(disp_array[:,0] - model_ssm))  #find simulation start time index
disp_dict['time']= disp_array[start_idx:,0] - model_ssm +1
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]-sfc_hgt))
disp_dict['meta']= 'time: seconds since model start run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'


#load pre-burn moisture profile
pre_burn_qv = np.genfromtxt(pre_moisture, skip_header=1, delimiter=',')

#get indecies of samples corresponding to model output times
tidx = [np.argmin(abs(disp_dict['time']-t)) for t in tsec] #times since start
dt = disp_dict['time'][1] - disp_dict['time'][0]

#construct KDtree from idealized grid
lat = nc_data.variables['XLAT'][0,:,:]
lon = nc_data.variables['XLONG'][0,:,:]
grid_coord = zip(lat.ravel(),lon.ravel())
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx][:,0:2])

#calculate H20 MR departure and extract point locations (assumes dry initial state)
qvapornan = np.copy(qinterp)
qvapornan[qvapornan<1e-30] = np.nan
mod_val, obs_val = np.empty((len(tsec))) * np.nan, np.empty((len(tsec))) * np.nan
obs_h = []
print('Finding nearest points....')
for nt in range(len(tsec)):
	print('...tstep: %s') %nt
	idxy,idxx = np.unravel_index(grid_id[nt],np.shape(lat))
	idxz = np.argmin(abs(lvl-disp_dict['lcn'][tidx][nt][2]))
	mod_val[nt] = qvapornan[nt,idxz,idxy,idxx]
	obs_val[nt] = disp_dict['CO2'][tidx[nt]]		#in percent by volume
	obs_h.append(disp_dict['lcn'][tidx][nt][2])

#================================EMISSIONS==================================
#extract and format emissions data
emis_dict = {}
emis_array = np.genfromtxt(emis_data, skip_header=1, usecols = [5,6,13,14], delimiter=',',skip_footer=emis_excl)

emis_dict['bkgd']= zip(emis_array[:,0] - model_ssm +1, emis_array[:,1] - model_ssm +1)
emis_dict['smoke'] = zip((emis_array[:,2] - model_ssm +1).astype(int), (emis_array[:,3] - model_ssm +1).astype(int))
emis_dict['meta']= 'bkgd: background slice start and end in sec from simulation start | \
					smoke: plume start and end in sec from simulation start | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation MSL'

#================================PLOTTING==================================
# #plot of CO2 slices
# plt.figure()
# plt.title('CO2 anomalies')
# plt.scatter(disp_dict['time'],disp_dict['CO2'])
# ax = plt.gca()
# for nSlice in range(len(emis_dict['smoke'])):
# 	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
# 	ax.fill_between(shade, 390,440, facecolor='gray', alpha=0.1, edgecolor='w')
# plt.ylim([390,440])
# plt.xlim([tsec[0],tsec[-1]])
# plt.show()


#plot of H2O slices
plt.figure()
plt.title('H2O')
plt.scatter(disp_dict['time'],disp_dict['H2O'])
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 0,1.5, facecolor='gray', alpha=0.1, edgecolor='w')
plt.ylim([0,1.5])
# plt.xlim([tsec[0],tsec[-1]])
plt.show()

#plot of model H20 slices overlayed with real emissions
plt.title('SIMULATED $Q_v$ ANOMALY ALONG FLIGHT PATH')
plt.scatter(disp_dict['time'][tidx], mod_val*1000000,c='red')
plt.plot(disp_dict['time'][tidx], mod_val*1000000,'b--')
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 0,20, facecolor='gray', alpha=0.3, edgecolor='w')
plt.ylim([0,15])
plt.xlim([0,tsec[-1]])
plt.xlabel('time [s]')
plt.ylabel('$H_{2}O$ mixing ratio anomaly [mg/kg]')
plt.tight_layout()
# plt.savefig(fig_dir + 'LES_Qv_flight_path.pdf')
plt.show()

#================================FLIGHT ANIMATION==================================
fig = plt.figure()
ax = p3.Axes3D(fig)

# create initial frame
point, = ax.plot([disp_dict['lcn'][0,0]],[disp_dict['lcn'][0,1]],[disp_dict['lcn'][0,2]], 'o')
ax.contourf(WLAT, WLONG, np.zeros(np.shape(WLAT)), alpha=0.3)
line, = ax.plot(disp_dict['lcn'][:,0], disp_dict['lcn'][:,1], disp_dict['lcn'][:,2], label='flight path', color='gray', alpha=0.3)
ax.legend()
ax.set_xlim([min(disp_dict['lcn'][:,0]), max(disp_dict['lcn'][:,0])])
ax.set_ylim([min(disp_dict['lcn'][:,1]), max(disp_dict['lcn'][:,1])])
ax.set_zlim([min(disp_dict['lcn'][:,2]), max(disp_dict['lcn'][:,2])])
time_text = ax.text(0.05,0.05,0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

#make a list of all times within plume from emissions
smoky = []
for item in emis_dict['smoke']:
	smoky.extend(np.arange(item[0],item[1]))

# move the point position at every frame
def update_point(n, disp_dict,smoky,point):
    point.set_data(np.array([disp_dict['lcn'][n,0],disp_dict['lcn'][n,1]]))
    point.set_3d_properties(disp_dict['lcn'][n,2], 'z')
    time_text.set_text('Time (sec) = %s' %(n*dt))
    if disp_dict['time'][n] in smoky:
    	point.set_color('r')
    else:
    	point.set_color('k')
    return point, time_text,

#plot the first 1500 frames (3000sec) - roughtly the length of the simulation
ani=animation.FuncAnimation(fig, update_point, 3000, fargs=(disp_dict,smoky,point), interval=15)
# ani.save('./test_ani.gif', writer='imagemagick',fps=120)
plt.show()


#================================VIRTICAL PROFILE==================================


g_s = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - garage_ssm[0])) 	#get index in disersion dict
g_f = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - garage_ssm[1]))

c_s = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - corskcrew_ssm[0]))
c_f = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - corskcrew_ssm[1]))


# s = np.argmax(disp_dict['lcn'][:,2]) 	#corkscrew starts at max height - get the time of max height after end of spinup
# f = s + 250 								#corkscrew lasts ~500sec (250x 2 sec timesteps)
# cleanf = np.argmin(disp_dict['lcn'][:,2])	#end of initial vertical profile (beginning of flight) - used as bg
# h20 = np.array(disp_dict['H2O'])
# h = np.array(disp_dict['lcn'])[:,2]
# bg_h2o = disp_dict['H2O'][:cleanf]
# bg_h = h[:cleanf]

# #calculate H20 MR departure at corkscrew point locations (NOT MATCHED IN TIME - LAST SLID ONLY)
# #H20 corkscrew profile from last frame
# dist_cs, grid_id_cs = gridTree.query(np.array(disp_dict['lcn'])[s:f,0:2])
# mod_val_cs = []
# print('Finding nearest points....')
# for nt in range(len(range(s,f))):
# 	print('...tstep: %s') %nt
# 	idxy,idxx = np.unravel_index(grid_id_cs[nt],np.shape(lat))
# 	idxz = np.argmin(abs(lvl-disp_dict['lcn'][nt,2]))
# 	mod_val_cs.append(qvapornan[-1,idxz,idxy,idxx]*1000000)
# 	# mv = disp_dict['H2O'][tidx[nt]]		#in percent by volume

# #define start and end of the corkscrew in sec from beginning of simulation
# plt.figure(figsize=(9,5))
# plt.subplot(1,2,1)
# plt.title('(a) $CO_2$ PROFILE FROM CORKSCREW')
# plt.scatter(disp_dict['CO2'][s:f],disp_dict['lcn'][s:f,2],color='black' )
# plt.plot(disp_dict['CO2'][s:f],disp_dict['lcn'][s:f,2],'k--' )
# plt.ylim([0,1700])
# plt.xlabel('$CO_2$ mixing ratio [ppmv]')
# plt.ylabel('height [m]')
# plt.subplot(1,2,2)
# plt.title('(b) $H_2O$ PROFILE FROM LES CORKSCREW')
# plt.scatter(mod_val_cs,disp_dict['lcn'][s:f,2],color='blue' )
# plt.plot(mod_val_cs,disp_dict['lcn'][s:f,2],'b--' )
# plt.ylim([0,1700])
# plt.xlabel('$H_2O$ mixing ratio [mg/kg]')
# plt.ylabel('height [m]')
# plt.tight_layout()
# plt.savefig(fig_dir + 'ProfilesCorkscrew.pdf')
# plt.show()


#H2O and CO2 profiles from garage flights
plt.figure(figsize=(9,5))
plt.subplot(1,2,1)
plt.title('(a) $CO_2$ PROFILE FROM CORKSCREW')
plt.scatter(obs_val[c_s:c_f],obs_h[c_s:c_f],color='black')
plt.plot(obs_val[c_s:c_f],obs_h[c_s:c_f],'k--' )
plt.ylabel('height [m]')
plt.ylim([0,1700])
plt.xlabel('$CO_2$ mixing ratio [ppmv]')
plt.subplot(1,2,2)
plt.title('(b) $H_2O$ PROFILE FROM LES CORKSCREW')
plt.scatter(mod_val[c_s:c_f]*1000000,obs_h[c_s:c_f],color='blue')
plt.plot(mod_val[c_s:c_f]*1000000,obs_h[c_s:c_f],'b--' )
plt.ylim([0,1700])
plt.xlabel('$H_2O$ mixing ratio [mg/kg]')
plt.ylabel('height [m]')
plt.tight_layout()
plt.savefig(fig_dir + 'ProfilesGarage.pdf')
plt.show()

#H2O and CO2 profiles from garage flights
plt.figure(figsize=(9,5))
plt.subplot(1,2,1)
plt.title('(a) $CO_2$ PROFILE FROM OBS GARAGE')
plt.scatter(obs_val[g_s:g_f],obs_h[g_s:g_f],color='black')
plt.plot(obs_val[g_s:g_f],obs_h[g_s:g_f],'k--' )
plt.ylabel('height [m]')
plt.ylim([0,1700])
plt.xlabel('$CO_2$ mixing ratio [ppmv]')
plt.subplot(1,2,2)
plt.title('(b) $H_2O$ PROFILE FROM LES GARAGE')
plt.scatter(mod_val[g_s:g_f]*1000000,obs_h[g_s:g_f],color='blue')
plt.plot(mod_val[g_s:g_f]*1000000,obs_h[g_s:g_f],'b--' )
plt.ylim([0,1700])
plt.xlabel('$H_2O$ mixing ratio [mg/kg]')
plt.ylabel('height [m]')
plt.tight_layout()
plt.savefig(fig_dir + 'ProfilesGarage.pdf')
plt.show()

#H2O profiles from corkscrew and earlier
plt.title('$H_2O$ PROFILES')
plt.plot(pre_burn_qv[:,1],pre_burn_qv[:,0],':', label='pre-burn $H_2O$ profile from sounding')
plt.plot(bg_h2o, bg_h,'g--',label='background $H_2O$ profile from flight ')
plt.plot(disp_dict['H2O'][s:f],disp_dict['lcn'][s:f,2],'r-.',label='in-plume $H_2O$ profile from flight' )
plt.ylim([0,1700])
plt.xlabel('$H_2O$ mixing ratio [%]')
plt.ylabel('height [m]')
plt.tight_layout()
plt.legend(loc='lower left')
plt.savefig(fig_dir + 'H2OProfiles.pdf')
plt.show()


#vertical column evoluation
column_evol = np.nansum(np.nansum(qinterp,2),2)
column_evol[column_evol<0] = np.nan 			#mask negataives
plt.pcolor(column_evol.T*1000,cmap=plt.cm.cubehelix_r)
cbar = plt.colorbar()
cbar.set_label('$H_{2}O$ mixing ratio anomaly [g/kg]')
ax = plt.gca()
ax.set_yticks(np.arange(0,numLvl,10))
ax.set_yticklabels(lvl[::10])
ax.set_xticks(np.arange(0,len(tsec),30))
ax.set_xticklabels((tsec[::30]/60.).astype(int))
plt.xlabel('time [min]')
plt.ylabel('height [m]')
plt.title('EVOLUTION OF SMOKE CONCENTRATION COLUMN')
plt.tight_layout()
plt.savefig(fig_dir + 'SmokeColumn.pdf')
plt.show()

# #create a horizontall averaged plume profile (last slide)
# plume_profile = np.nansum(qinterp[-1,:,:,:],2)
# plt.pcolor(plume_profile * 1000,cmap=plt.cm.PuBu)
# cbar = plt.colorbar()
# cbar.set_label('$H_{2}O$ mixing ratio anomaly [g/kg]')
# ax = plt.gca()
# ax.set_yticks(np.arange(0,numLvl,5))
# ax.set_yticklabels(lvl[::5])
# plt.xlabel('grid #')
# plt.ylabel('height [m]')
# plt.title('TOTAL $Q_v$ ANOMALY (HORIZONTALLY AVERAGED) at $t_{end}$')
# plt.tight_layout()
# plt.savefig(fig_dir + 'AvePlume_Tend.pdf')
# plt.show()


# #=============================ANIMATION OF CW ave plume==================================

# print('WARNING: Slow routine: rotating the array to be alighned with mean wind')
# qrot = rotate(qinterp[:,:,:,:], 40, axes=(2, 3), reshape=True, mode='constant', cval=np.nan)
# qrot[qrot<1e-30] = np.nan
# print('WARNING: Slow routine: creating cross-wind averages')
# cw_sum = np.nansum(qrot,3)*1000 #converting to mg

# fig = plt.figure()
# ax = plt.gca()
# # create initial frame
# # point, = ax.plot([disp_dict['lcn'][0,0]],[disp_dict['lcn'][0,1]],[disp_dict['lcn'][0,2]], 'o')
# cntr = plt.contourf(cw_sum[0,:,:],cmap=plt.cm.PuBu,levels=np.arange(0,0.9,0.1))
# cbar = plt.colorbar()
# cbar.set_label('total $H_{2}O$ mixing ratio anomaly [g/kg]')
# plt.xlabel('grid #')
# plt.ylabel('height [m]')
# ax.set_yticks(np.arange(0,numLvl,5))
# ax.set_yticklabels(lvl[::5])
# plt.title('EVOLUTION OF TOTAL CROSS-WIND $Q_v$ ANOMALY')

# # plt.clim([0,2])

# # line, = ax.plot(disp_dict['lcn'][:,0], disp_dict['lcn'][:,1], disp_dict['lcn'][:,2], label='flight path', color='gray', alpha=0.3)
# # ax.legend()
# # ax.set_xlim([min(disp_dict['lcn'][:,0]), max(disp_dict['lcn'][:,0])])
# # ax.set_ylim([min(disp_dict['lcn'][:,1]), max(disp_dict['lcn'][:,1])])
# # ax.set_zlim([min(disp_dict['lcn'][:,2]), max(disp_dict['lcn'][:,2])])
# # ax.colorbar()
# # time_text = ax.text(0.05,0.05,0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

# # move the point position at every frame
# def update_plot(n, cw_ave,cntr):
#     cntr = plt.contourf(cw_sum[n,:,:],cmap=plt.cm.PuBu,levels=np.arange(0,0.9,0.1))
#     # time_text.set_text('Time (sec) = %s' %(n*dt))
#     return cntr,

# #plot the first 1500 frames (3000sec) - roughtly the length of the simulation
# ani=animation.FuncAnimation(fig, update_plot, 199, fargs=(cw_sum,cntr), interval=1)
# ani.save(fig_dir + 'CW_total.gif', writer='imagemagick',fps=120)
# plt.show()



#NEED TO ADD COLORBAR AND CONSTANT LIMITS
