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
import pyproj
import os.path
import pickle
import mpl_toolkits.mplot3d.axes3d as p3
import mpl_toolkits.mplot3d as a3
from matplotlib import animation
import pandas as pd
import matplotlib.dates as mdates
import datetime as dt
import imp
import wrf

#====================INPUT===================
#all common variables are stored separately
import rxcadreDRY as rx
imp.reload(rx)		        #force load each time

# emis_excl = 0 							#number of samples to excluded (from the END!)
rskcrew_ssm= [47183,47521]			#start and end of rskcrew maneuver in ssm
bg_rk_ssm = [43975,44371] 			#start and end time of pre-burn corkscrew for background
garage_ssm = [45000,47200] 				#start and end time of garage profile

animations = 0
#=================end of input===============
print('ANALYSIS OF VERTICAL PLUME RISE AND DIPSERSION')

print('.....extracting NetCDF data from %s ' %rx.wrfdata)
ncdata = netcdf.netcdf_file(rx.wrfdata, mode ='r')

#load georeferencing data for the same run
print('.....importing wrf coordinates from  %s ' %rx.geo_path)
wrfgeo = np.load(rx.geo_path, allow_pickle=True).item()

#get variables
ncdict = wrf.extract_vars(ncdata, None, ('PHB','PH','XTIME','tr17_1','tr17_2'))

#get height and destagger vars
print('.....destaggering height data')
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)

#get domain dimensions
nT,nZ,nY,nX = np.shape(z)

#open/generate basemap
if os.path.isfile(rx.basemap_path):
    bm = pickle.load(open(rx.basemap_path,'rb'))   # load here the above pickle
    print('Domain basemap found at: %s' %rx.basemap_path)
else:
    print('WARNING: no existing basemaps found: configuring a new basemap')
    bm = basemap.Basemap(llcrnrlon=wrfgeo['WLONG'][0,0], llcrnrlat=wrfgeo['WLAT'][0,0],\
                            urcrnrlon=wrfgeo['WLONG'][-1,-1], urcrnrlat=wrfgeo['WLAT'][-1,-1], resolution='f', epsg=4326)
    pickle.dump(bm,open(rx.basemap_path,'wb'),-1)  	# pickle the new map for later
    print('.....New basemap instance saved as: %s' %rx.basemap_path)

# # Sanity check: import shape file
# polygons = bm.readshapefile(rx.bounds_shape,name='fire_bounds',drawbounds=True)
# fireim = ncdata.variables['GRNHFX'][200,:,:]
# bm.contourf(wrfgeo['WLONG'],wrfgeo['WLAT'],fireim)
# plt.show()

#extract model time info
tsec = ncdict['XTIME'][:] * 60. 		#get time in seconds since spinup start
model_ssm = int(rx.runstart[0:2])*3600 + int(rx.runstart[3:5])*60

#==========================VERTICAL INTERPOLATION============================
numLvl = len(rx.lvl)

#open/generate vertically interpolated data
if os.path.isfile(rx.interp_path):
    qinterp = np.load(rx.interp_path)   # load here the above pickle
    print('Interpolated data found at: %s' %rx.interp_path)
else:
    print('WARNING: no interpolated vapour data found - generating: SLOW ROUTINE!')
    qvcopy = np.copy(ncdata.variables['QVAPOR'][:,:,:,:])
    qinterp = np.empty((nT,numLvl,nY,nX)) * np.nan
    for t in range(nT):
        print('.... tsetp = %s/%s' %(t,nT))
        for y in range(nY):
            for x in range(nX):
                tempz = z[t,:,y,x]
                interpz = (tempz[:-1]+tempz[1:])/2.
                f = interpolate.interp1d(interpz,qvcopy[t,:,y,x],fill_value="extrapolate")
                qinterp[t,:,y,x] = f(rx.lvl)
                print('.... masking cells with vals < 1e-30')
                qinterp[qinterp<1e-30] = np.nan
                np.save(rx.interp_path, qinterp)
                print('Interpolated data saved as: %s' %rx.interp_path)

#================================DISPERSION==================================
#extract and format dispersion data
print('Importing dispersion data from %s' %rx.disp_data)
disp_dict = {}
disp_array = np.genfromtxt(rx.disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')

start_idx = np.argmin(abs(disp_array[:,0] - model_ssm))  #find simulation start time index
disp_dict['time']= disp_array[start_idx:,0] - model_ssm + 1
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]-rx.sfc_hgt))
disp_dict['meta']= 'time: seconds since restart run | \
    CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
    CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
    CH4: Mixing ratio of methane in units of ppmv in dry air. | \
    H2O: Mixing ratio of water vapor in percent by volume. | \
    lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'


#load pre-burn moisture profile
pre_burn_qv = np.genfromtxt(rx.pre_moisture, skip_header=1, delimiter=',')

#get indecies of samples corresponding to model output times
tidx = [np.argmin(abs(disp_dict['time']-t)) for t in tsec] #times since start
delt = disp_dict['time'][1] - disp_dict['time'][0]


#convert all timestamps to datetime objects
basetime = dt.datetime(year=2012,month=11,day=10)
timestamp = [basetime + dt.timedelta(hours = 10, seconds = i) for i in disp_dict['time']]


#construct KDtree from idealized grid
lat = ncdata.variables['XLAT'][0,:,:]
lon = ncdata.variables['XLONG'][0,:,:]
grid_coord = zip(lat.ravel(),lon.ravel())
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx][:,0:2])

#calculate H20 MR departure and extract point locations (assumes dry initial state)

# qvapornan = np.copy(qinterp) 			#mask empty cells
# qvapornan[qvapornan<1e-30] = np.nan
mod_val, obs_val = np.empty((len(tsec))) * np.nan, np.empty((len(tsec))) * np.nan
obs_h = []
obs_lon, obs_lat = [],[]
print('Finding nearest points....')
for nt in range(len(tsec)):
    print('...tstep: %s') %nt
    idxy,idxx = np.unravel_index(grid_id[nt],np.shape(lat))
    idxz = np.argmin(abs(rx.lvl-disp_dict['lcn'][tidx][nt][2]))
    mod_val[nt] = qinterp[nt,idxz,idxy,idxx]
    obs_val[nt] = disp_dict['CO2'][tidx[nt]]		#in percent by volume
    obs_h.append(disp_dict['lcn'][tidx][nt][2])
    obs_lon.append(disp_dict['lcn'][tidx][nt][1])
    obs_lat.append(disp_dict['lcn'][tidx][nt][0])
#================================EMISSIONS==================================
print('Importing emissions data from %s' %rx.emis_data)
#extract and format emissions data
emis_dict = {}
emis_array = np.genfromtxt(rx.emis_data, skip_header=1, usecols = [5,6,13,14], delimiter=',',skip_footer=emis_excl)

emis_dict['bkgd']= zip(emis_array[:,0] - model_ssm +1, emis_array[:,1] - model_ssm +1)
emis_dict['smoke'] = zip((emis_array[:,2] - model_ssm +1).astype(int), (emis_array[:,3] - model_ssm +1).astype(int))
emis_dict['meta']= 'bkgd: background slice start and end in sec from simulation start | \
    smoke: plume start and end in sec from simulation start | \
    lcn: (lat, lon, elevation) - coords in WGS84, elevation MSL'

#================================PLOTTING==================================

#plot of model H20 slices overlayed with real emissions
plt.title('SIMULATED $Q_v$ ANOMALY ALONG FLIGHT PATH')
plt.plot(np.array(timestamp)[tidx], mod_val*1000000,'ro')
plt.plot(np.array(timestamp)[tidx], mod_val*1000000,'b--')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
for nSlice in range(len(emis_dict['smoke'])):
    smoke_start = basetime + dt.timedelta(hours = 10, seconds = emis_dict['smoke'][nSlice][0])
    smoke_end =  basetime + dt.timedelta(hours = 10, seconds = emis_dict['smoke'][nSlice][1])
    ax.fill_between(pd.date_range(smoke_start,smoke_end,freq='S'),0,20,facecolor='gray', alpha=0.3)
plt.ylim([0,15])
plt.gcf().autofmt_xdate()
plt.xlim([dt.datetime(year=2012,month=11,day=10,hour=12,minute=27),dt.datetime(year=2012,month=11,day=10,hour=13,minute=4)])   #hardcoded
plt.xlabel('time (CST)')
plt.ylabel('$H_{2}O$ mixing ratio anomaly [mg/kg]')
plt.tight_layout()
plt.savefig(rx.fig_dir + 'LES_Qv_flight_path.pdf')
plt.show()
#
#================================FLIGHT ANIMATION==================================
# fig = plt.figure()
# ax = p3.Axes3D(fig)
# ax = plt.gca()
# # create initial frame
# point, = ax.plot([disp_dict['lcn'][0,0]],[disp_dict['lcn'][0,1]],[disp_dict['lcn'][0,2]], 'o')
# ax.contourf(WLAT, WLONG, np.zeros(np.shape(WLAT)), alpha=0.3)
# line, = ax.plot(disp_dict['lcn'][:,0], disp_dict['lcn'][:,1], disp_dict['lcn'][:,2], label='flight path', color='gray', alpha=0.3)
# ax.legend()
# ax.set_xlim([min(disp_dict['lcn'][:,0]), max(disp_dict['lcn'][:,0])])
# ax.set_ylim([min(disp_dict['lcn'][:,1]), max(disp_dict['lcn'][:,1])])
# ax.set_zlim([min(disp_dict['lcn'][:,2]), max(disp_dict['lcn'][:,2])])
# time_text = ax.text(0.05,0.05,0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
#
# #make a list of all times within plume from emissions
# smoky = []
# for item in emis_dict['smoke']:
# 	smoky.extend(np.arange(item[0],item[1]))
#
# # move the point position at every frame
# def update_point(n, disp_dict,smoky,point):
#     point.set_data(np.array([disp_dict['lcn'][n,0],disp_dict['lcn'][n,1]]))
#     point.set_3d_properties(disp_dict['lcn'][n,2], 'z')
#     time_text.set_text('Time (sec) = %s' %(n*delt))
#     if disp_dict['time'][n] in smoky:
#     	point.set_color('r')
#     else:
#     	point.set_color('k')
#     return point, time_text,
#
# #plot the first 1500 frames (3000sec) - roughtly the length of the simulation
# ani=animation.FuncAnimation(fig, update_point, 1500, fargs=(disp_dict,smoky,point), interval=15)
# # ani.save('./test_ani.gif', writer='imagemagick',fps=120)
# plt.show()
# plt.close()

if animations:
    print('Animating top view of flight....')
    fig = plt.figure()

    # create initial frame
    smokeim = np.nansum(qinterp[0,:,:,:],0) * 1000
    im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower')
    scat = bm.scatter(disp_dict['lcn'][0,1],disp_dict['lcn'][0,0],40,marker='o')

    #make a list of all times within plume from emissions
    smoky = []
    for item in emis_dict['smoke']:
        smoky.extend(np.arange(item[0],item[1]))

    # move the point position at every frame
    def update_point(n, disp_dict,smoky,scat,im):
        m = tidx[n]
        # i = int(np.floor(n/5.))
        del im
        smokeim = np.nansum(qinterp[n,:,:,:],0) * 1000
        im = bm.imshow(smokeim, cmap = plt.cm.bone_r)
        # im.set_data(smokeim)
        scat.set_offsets(np.c_[disp_dict['lcn'][m,1],disp_dict['lcn'][m,0]])
        # plt.gca().set_title('')

        if disp_dict['time'][m] in smoky:
            scat.set_color('r')
        else:
            scat.set_color('k')

    #plot the first 1500 frames (3000sec) - roughtly the length of the simulation
    # ani=animation.FuncAnimation(fig, update_point, 1349, fargs=(disp_dict,smoky,scat,im), interval=10)
    ani=animation.FuncAnimation(fig, update_point, 269, fargs=(disp_dict,smoky,scat,im),interval = 100, repeat=0)

    ani.save(rx.fig_dir + 'FlightTopView.gif', writer='imagemagick',fps=120)
    # plt.show()
    plt.close()


#================================VIRTICAL PROFILE==================================


g_s = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - garage_ssm[0])) 	#get index in dispersion dict
g_f = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - garage_ssm[1]))

c_s = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - corskcrew_ssm[0]))
c_f = np.argmin(abs(disp_dict['time'][tidx] + model_ssm - corskcrew_ssm[1]))
dict_c_s = np.argmin(abs(disp_dict['time'] + model_ssm - corskcrew_ssm[0]))
dict_c_f = np.argmin(abs(disp_dict['time'] + model_ssm - corskcrew_ssm[1]))
dict_bg_s = np.argmin(abs(disp_dict['time'] + model_ssm - bg_cork_ssm[0]))
dict_bg_f = np.argmin(abs(disp_dict['time'] + model_ssm - bg_cork_ssm[1]))
# dict_bg_s = np.argmin(abs(disp_dict['time'] + model_ssm - garage_ssm[0]))
# dict_bg_f = np.argmin(abs(disp_dict['time'] + model_ssm - garage_ssm[1]))

# T0 = ncdata.variables['T'][0,:,:,:]
# Tprofile = T0[:,0,0]
# plt.plot(Tprofile, z[0,:-1,0,0])
# plt.show()




#top view of smoke for corkscrew with average obs location
cs_lon =  np.mean(disp_dict['lcn'][dict_c_s:dict_c_f,1])
cs_lat = np.mean(disp_dict['lcn'][dict_c_s:dict_c_f,0])
smokeim = np.nansum(qinterp[258,:,:,:],0) * 1000
im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower')
bm.scatter(cs_lon,cs_lat,40,marker='*',color='r')
plt.colorbar(im, label='total column $H_{2}O$ mixing ratio anomaly [mg/kg]')
plt.tight_layout()
plt.savefig(rx.fig_dir + 'CSLocation.pdf')
plt.show()


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
plt.savefig(rx.fig_dir + 'ProfilesCorkscrew.pdf')
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
plt.savefig(rx.fig_dir + 'ProfilesGarage.pdf')
plt.show()
#
# #H2O profiles from corkscrew and earlier
# plt.title('$H_2O$ PROFILES')
# plt.plot(pre_burn_qv[:,1],pre_burn_qv[:,0],':', label='pre-burn $H_2O$ profile from sounding')
# plt.plot(disp_dict['H2O'][dict_bg_s:dict_bg_f],disp_dict['lcn'][dict_bg_s:dict_bg_f,2],'g--',label='pre-burn $H_2O$ profile from flight ')
# plt.plot(disp_dict['H2O'][dict_c_s:dict_c_f],disp_dict['lcn'][dict_c_s:dict_c_f,2],'r-.',label='in-plume $H_2O$ profile from flight' )
# plt.ylim([0,1700])
# plt.xlabel('$H_2O$ mixing ratio [%]')
# plt.ylabel('height [m]')
# plt.tight_layout()
# plt.legend(loc='lower left')
# plt.savefig(rx.fig_dir + 'H2OProfiles.pdf')
# plt.show()


#vertical column evoluation
column_evol = np.nansum(np.nansum(qinterp,2),2)
column_evol[column_evol<0] = np.nan 			#mask negataives
plt.pcolor(column_evol.T*1000,cmap=plt.cm.cubehelix_r)
cbar = plt.colorbar()
cbar.set_label('$H_{2}O$ mixing ratio anomaly [g/kg]')
ax = plt.gca()
ax.set_yticks(np.arange(0,numLvl,10))
ax.set_yticklabels(rx.lvl[::10])
# ax.xaxis_date()
ax.set_xticks(np.arange(0,len(tsec),50))
ax.set_xticklabels([np.array(timestamp)[tidx][i].strftime('%H:%M') for i in np.arange(0,len(tsec),50)])
plt.xlabel('time [CST]')
plt.ylabel('height [m]')
plt.title('EVOLUTION OF SMOKE CONCENTRATION COLUMN')
plt.tight_layout()
plt.savefig(rx.fig_dir + 'SmokeColumn.pdf')
plt.show()

#top view of smoke for corkscrew with average obs location
cs_lon =  np.mean(disp_dict['lcn'][dict_c_s:dict_c_f,1])
cs_lat = np.mean(disp_dict['lcn'][dict_c_s:dict_c_f,0])
plt.title('TOP VIEW CORKSCREW')
smokeim = np.nansum(qinterp[258,:,:,:],0) * 1000
im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower')
bm.scatter(cs_lon,cs_lat,40,marker='*',color='r')
plt.colorbar(im, label='total column $H_{2}O$ mixing ratio anomaly [mg/kg]')
plt.tight_layout()
plt.savefig(rx.fig_dir + 'CSLocation.pdf')
plt.show()


#=============================ANIMATION OF CW ave plume==================================
if animations:
    print('WARNING: Slow routine: rotating the array to be alighned with mean wind')
    qzero = qinterp*1000
    qzero[np.isnan(qzero)] = 0

    cw_sum = []
    for nTime in range(nT):
        print(nTime)
        qrot = rotate(qzero[nTime,:,:,:], 39, axes=(1, 2), reshape=True, mode='constant', cval=np.nan)
        totq = np.nansum(qrot,2)
        cw_sum.append(totq)

    cw_sum = np.array(cw_sum)
    fig = plt.figure()
    ax = plt.gca()
    # create initial frame
    cntr = plt.contourf(cw_sum[0,:,:],cmap=plt.cm.PuBu,levels=np.arange(0,0.9,0.1))
    cbar = plt.colorbar()
    cbar.set_label('total $H_{2}O$ mixing ratio anomaly [mg/kg]')
    plt.xlabel('distance [km]')
    plt.ylabel('height [m]')
    ax.set_yticks(np.arange(0,numLvl,5))
    ax.set_yticklabels(rx.lvl[::5])
    ax.set_xticks(np.arange(0,450,50))
    ax.set_xticklabels((np.arange(0,450,50)*0.040).astype(int))
    plt.title('EVOLUTION OF TOTAL CROSS-WIND $Q_v$ ANOMALY')

    # update the frame
    def update_plot(n,cw_sum,cntr):
        del cntr
        cntr = plt.contourf(cw_sum[n,:,:],cmap=plt.cm.PuBu,levels=np.arange(0,0.9,0.1))
        return cntr,

    #plot all frames
    ani=animation.FuncAnimation(fig, update_plot, 270, fargs=(cw_sum,cntr), interval=1)
    ani.save(rx.fig_dir + 'Plume_CS_animation.gif', writer='imagemagick',fps=120)
    # plt.show()
