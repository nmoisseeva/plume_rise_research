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
import sys

#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
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
ncdict['CO2'] = ncdict.pop('tr17_1')
ncdict['CO'] = ncdict.pop('tr17_2')
tracers = ['CO','CO2']

#get height and destagger
if os.path.isfile(rx.z_path):
    print('.....loading destaggered height array from %s' %rx.z_path)
    z = np.load(rx.z_path)
else:
    print('.....destaggering height data')
    zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
    z = wrf.destagger(zstag,1)
    np.save(rx.z_path, z)

#get domain dimensions
nT,nZ,nY,nX = np.shape(z)

#open/generate basemap
if os.path.isfile(rx.basemap_path):
    bm = pickle.load(open(rx.basemap_path,'rb'))   # load here the above pickle
    print('.....domain basemap found at: %s' %rx.basemap_path)
else:
    print('WARNING: no existing basemaps found: configuring a new basemap')
    bm = basemap.Basemap(llcrnrlon=wrfgeo['WLONG'][0,0], llcrnrlat=wrfgeo['WLAT'][0,0],\
                            urcrnrlon=wrfgeo['WLONG'][-1,-1], urcrnrlat=wrfgeo['WLAT'][-1,-1], resolution='f', epsg=4326)
    pickle.dump(bm,open(rx.basemap_path,'wb'),-1)  	# pickle the new map for later
    print('.....new basemap instance saved as: %s' %rx.basemap_path)

# # Sanity check: import shape file
# polygons = bm.readshapefile(rx.bounds_shape,name='fire_bounds',drawbounds=True)
# fireim = ncdata.variables['GRNHFX'][200,:,:]
# bm.contourf(wrfgeo['WLONG'],wrfgeo['WLAT'],fireim)
# plt.show()

#extract model time info
model_ssm = int(rx.runstart[0:2])*3600 + int(rx.runstart[3:5])*60

#==========================VERTICAL INTERPOLATION============================
numLvl = len(rx.lvl)

#open/generate vertically interpolated data
for tracer in tracers:
    tracer_path = '%s%s.npy' %(rx.interp_path, tracer)
    if os.path.isfile(tracer_path):
        print('.....interpolated %s data found at: %s' %(tracer,tracer_path))
    else:
        print('WARNING: no interpolated smoke data found - generating: SLOW ROUTINE!')
        tracerinterp = np.empty((nT,numLvl,nY,nX))
        print('.....processing tracer: %s' %tracer)
        sys.stdout.write("[%s]" % (" " * nT))
        sys.stdout.flush()
        sys.stdout.write("\b" * (nT+1)) # setup progress bar for loop
        for t in range(nT):
            for y in range(nY):
                for x in range(nX):
                    f = interpolate.interp1d(z[t,:,y,x],ncdict[tracer][t,:,y,x],fill_value="extrapolate")
                    interpvector = f(rx.lvl)
                    # interpvector[interpvector<1e-30] = np.nan           #mask tiny and negative values
                    tracerinterp[t,:,y,x] = interpvector
            sys.stdout.write("-")
            sys.stdout.flush()
        sys.stdout.write("]\n") # end progress bar
        np.save(tracer_path, tracerinterp)
        print('Interpolated data saved as: %s' %tracer_path)

#================================DISPERSION==================================
basetime = dt.datetime(year=2012,month=11,day=10)

#extract and format dispersion data
print('.....importing dispersion data from %s' %rx.disp_data)
disp_dict = {}
disp_array = np.genfromtxt(rx.disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')
disp_dict['time'] = np.array([basetime + dt.timedelta(seconds = int(i)) for i in disp_array[:,0]])
disp_dict['CO'] = disp_array[:,1]
disp_dict['CO2'] = disp_array[:,2]
disp_dict['CH4'] = disp_array[:,3]
disp_dict['H2O'] = disp_array[:,4]
disp_dict['lcn'] = np.array(list(zip(disp_array[:,5],disp_array[:,6],disp_array[:,7]-rx.sfc_hgt)))
disp_dict['meta']= 'time: min since restart run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'
# #load pre-burn moisture profile
# pre_burn_qv = np.genfromtxt(rx.pre_moisture, skip_header=1, delimiter=',')

#define timing of events
model_datetime = [basetime + dt.timedelta(hours=int(rx.runstart[:2]),minutes=t) for t in ncdict['XTIME']]
cs_start, cs_end = dt.datetime.combine(basetime,dt.datetime.strptime(rx.corkscrew[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(rx.corkscrew[1], '%H:%M:%S').time())
gg_start, gg_end = dt.datetime.combine(basetime,dt.datetime.strptime(rx.garage[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(rx.garage[1], '%H:%M:%S').time())

#get indecies of samples corresponding to model output times
tidx = [np.argmin(abs(disp_dict['time'] - t)) for t in model_datetime] #times since start

#construct KDtree from model grid
grid_coord = list(zip(wrfgeo['WLAT'].ravel(),wrfgeo['WLONG'].ravel()))
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx][:,0:2])

#extract model values for CO2 and CO along flight path and save with corresponding observed samples
model_dict ={'CO': np.empty((nT)), 'CO2': np.empty((nT)), 'Z':np.empty((nT))}
obs_dict = {'CO': np.empty((nT)), 'CO2': np.empty((nT)), 'Z':np.empty((nT))}
molar_mass = {'CO': 28.0101, 'CO2': 44.0095, 'air': 28.9644}
for tracer in tracers:
    tracer_path = '%s%s.npy' %(rx.interp_path, tracer)
    print('.....extracting data along flight path from %s' %tracer_path)
    tracerinterp = np.load(tracer_path)
    for nt in range(nT):
        idxy,idxx = np.unravel_index(grid_id[nt],np.shape(wrfgeo['XLONG']))
        idxz = np.argmin(abs(rx.lvl - disp_dict['lcn'][tidx][nt][2]))
        tracer_val_mukg = tracerinterp[nt,idxz,idxy,idxx]
        model_dict[tracer][nt] = 1e-3 * tracer_val_mukg * molar_mass['air'] / molar_mass[tracer]    #convert to ppmv, assuming model is in mug/kg
        obs_dict[tracer][nt] = disp_dict[tracer][tidx[nt]]
        obs_dict['Z'][nt] = disp_dict['lcn'][tidx][nt][2]
        model_dict['Z'][nt] = rx.lvl[idxz]

#================================EMISSIONS==================================
print('.....importing emissions data from %s' %rx.emis_data)
#extract and format emissions data
emis_dict = {}
emis_array = np.genfromtxt(rx.emis_data, skip_header=1, usecols = [4, 7,8,13,14], delimiter=',')
emis_dict['smoke_start'] = np.array([basetime + dt.timedelta(seconds = int(i)) for i in emis_array[:,3]])
emis_dict['smoke_end'] = np.array([basetime + dt.timedelta(seconds = int(i)) for i in emis_array[:,4]])
emis_dict['bkgdCO2'] = emis_array[:,1]
emis_dict['bkgdCO'] = emis_array[:,2]
emis_dict['elevation'] = emis_array[:,0]
emis_dict['meta']= 'smoke start/end: plume entry and exit points; background concentrations'

#interpolated background values for all levels
bginterpCO = interpolate.interp1d(emis_dict['elevation'],emis_dict['bkgdCO'],fill_value="extrapolate")
bginterpCO2 = interpolate.interp1d(emis_dict['elevation'],emis_dict['bkgdCO2'],fill_value="extrapolate")

#================================PLOTTING==================================

#plot of model tracers overlayed with real emissions
plt.title('CO ALONG FLIGHT PATH')
# plt.plot(disp_dict['time'][tidx], model_dict['CO'], 'r')
plt.plot(disp_dict['time'][tidx], model_dict['CO'], 'r--', label='WRF-SFIRE CO concentrations')
plt.plot(disp_dict['time'][tidx], obs_dict['CO'] - bginterpCO(obs_dict['Z']), 'k', label='observed CO concentrations - background')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
for nSlice in range(len(emis_dict['smoke_start'])):
    ax.fill_between(pd.date_range(emis_dict['smoke_start'][nSlice],emis_dict['smoke_end'][nSlice],freq='S'),0,1.7,facecolor='gray', alpha=0.3)
# plt.ylim([0,15])
plt.gcf().autofmt_xdate()
plt.xlim([model_datetime[36],model_datetime[-42]])
plt.ylim([0,1.7])
plt.xlabel('time (CST)')
plt.ylabel('CO concentration [ppmv]')
plt.legend()
plt.tight_layout()
plt.savefig(rx.fig_dir + 'LES_CO_flight_path.pdf')
plt.show()
#
#================================FLIGHT ANIMATION==================================
# if animations:
#     print('Animating top view of flight....')
#     fig = plt.figure()
#
#     # create initial frame
#     smokeim = np.nansum(qinterp[0,:,:,:],0) * 1000
#     im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower')
#     scat = bm.scatter(disp_dict['lcn'][0,1],disp_dict['lcn'][0,0],40,marker='o')
#
#     #make a list of all times within plume from emissions
#     smoky = []
#     for item in emis_dict['smoke']:
#         smoky.extend(np.arange(item[0],item[1]))
#
#     # move the point position at every frame
#     def update_point(n, disp_dict,smoky,scat,im):
#         m = tidx[n]
#         # i = int(np.floor(n/5.))
#         del im
#         smokeim = np.nansum(qinterp[n,:,:,:],0) * 1000
#         im = bm.imshow(smokeim, cmap = plt.cm.bone_r)
#         # im.set_data(smokeim)
#         scat.set_offsets(np.c_[disp_dict['lcn'][m,1],disp_dict['lcn'][m,0]])
#         # plt.gca().set_title('')
#
#         if disp_dict['time'][m] in smoky:
#             scat.set_color('r')
#         else:
#             scat.set_color('k')
#
#     #plot the first 1500 frames (3000sec) - roughtly the length of the simulation
#     # ani=animation.FuncAnimation(fig, update_point, 1349, fargs=(disp_dict,smoky,scat,im), interval=10)
#     ani=animation.FuncAnimation(fig, update_point, 269, fargs=(disp_dict,smoky,scat,im),interval = 100, repeat=0)
#
#     ani.save(rx.fig_dir + 'FlightTopView.gif', writer='imagemagick',fps=120)
#     # plt.show()
#     plt.close()

#================================VIRTICAL PROFILE==================================
#get indecies of maneuver times
ics = np.argmin(abs(disp_dict['time'][tidx] - cs_start))
fcs = np.argmin(abs(disp_dict['time'][tidx] - cs_end))
igg = np.argmin(abs(disp_dict['time'][tidx] - gg_start))
fgg = np.argmin(abs(disp_dict['time'][tidx] - gg_end))

#top view of smoke for corkscrew with average obs location

#HARDCODED:
rotCS_lcn = [30.56683,-86.75389]    #calculated by rotating 30
ROT_cs_dist, ROT_cs_grid_id = gridTree.query(np.array(rotCS_lcn))
ROTidxy,ROTidxx = np.unravel_index(ROT_cs_grid_id,np.shape(wrfgeo['XLONG']))
ROT_profile_mukg = tracerinterp[258,:,ROTidxy,ROTidxx]
ROT_profile = 1e-3 * ROT_profile_mukg * molar_mass['air'] / molar_mass['CO2']    #convert to ppmv, assuming model is in mug/kg
#------end hardcoding

cs_lon =  np.mean(disp_dict['lcn'][tidx][ics:fcs,1])
cs_lat = np.mean(disp_dict['lcn'][tidx][ics:fcs,0])
tot_column_mukg = np.nansum(ncdict['CO'][258,:,:,:],0)
smokeim = 1e-3 * tot_column_mukg * molar_mass['air']/molar_mass['CO']   #convert to ppmv
im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower',vmin=0,vmax=50)
bm.scatter(cs_lon,cs_lat,40,marker='*',color='r',label='average location of corkscrew profile')
bm.scatter(rotCS_lcn[1],rotCS_lcn[0],20,marker='o',color='k',label='wind-corrected profile location')
plt.colorbar(im, label='total column CO [ppmv]')
plt.legend()
plt.tight_layout()
plt.savefig(rx.fig_dir + 'CSLocation.pdf')
plt.show()
plt.close()


#CO2 profiles from garage flights
plt.figure(figsize=(10,6))
plt.subplot(1,2,1)
plt.title('(a) $CO_2$ PROFILE FROM GARAGE PROFILE')
plt.plot(obs_dict['CO2'][igg:fgg]- bginterpCO2(obs_dict['Z'][igg:fgg]),obs_dict['Z'][igg:fgg],'k', label='observed CO$_2$ - background' )
plt.plot(model_dict['CO2'][igg:fgg] ,model_dict['Z'][igg:fgg],'r--', label='WRF-SFIRE CO$_2$')
plt.ylabel('height [m]')
plt.xlabel('$CO_2$ mixing ratio [ppmv]')
plt.ylim([0,1700])

plt.subplot(1,2,2)
plt.title('(b) $CO_2$ PROFILE FROM CORKSCREW PROFILE')
plt.plot(obs_dict['CO2'][ics:fcs] - bginterpCO2(obs_dict['Z'][ics:fcs]),obs_dict['Z'][ics:fcs],'k',label='observed CO$_2$ - background' )
plt.plot(model_dict['CO2'][ics:fcs], model_dict['Z'][ics:fcs],'r--',label='WRF-SFIRE CO$_2$')
plt.plot(ROT_profile,rx.lvl,'r:', label='rotated WRF-SFIRE CO$_2$')
plt.ylabel('height [m]')
plt.xlabel('$CO_2$ mixing ratio [ppmv]')
plt.ylim([0,1700])
plt.legend(loc='lower right')

plt.tight_layout()
plt.savefig(rx.fig_dir + 'CO2Profiles.pdf')
plt.show()
plt.close()


#vertical column evoluation
column_mukg = np.sum(np.sum(tracerinterp,2),2)
column_ppmv = 1e-3 * column_mukg * molar_mass['air']/molar_mass['CO2']
plt.pcolor(column_ppmv.T,cmap=plt.cm.cubehelix_r)
plt.colorbar(label='total domain CO$_2$ anomaly [ppmv]')
ax = plt.gca()
plt.ylabel('height [m]')
ax.set_yticks(np.arange(0,numLvl,10))
ax.set_yticklabels(rx.lvl[::10])
plt.xlabel('time [CST]')
ax.set_xticks(np.arange(0,nT,40))
ax.set_xticklabels([np.array(model_datetime)[i].strftime('%H:%M') for i in np.arange(0,nT,40)])
plt.title('EVOLUTION OF SMOKE CONCENTRATION COLUMN')
plt.tight_layout()
plt.savefig(rx.fig_dir + 'SmokeColumn.pdf')
plt.show()
plt.close()

#=============================ANIMATION OF CW ave plume==================================
if animations:
    print('WARNING: Slow routine: rotating the array to be alighned with mean wind')

    cw_sum = []
    for nTime in range(nT):
        print(nTime)
        CO2rot = rotate(tracerinterp[nTime,:,:,:], 39, axes=(1, 2), reshape=True, mode='constant')
        CO2tot = np.sum(CO2rot,2)
        cw_sum.append(CO2tot)

    cw_sum = np.array(cw_sum)
    fig = plt.figure()
    ax = plt.gca()
    # create initial frame
    cntr = plt.contourf(cw_sum[0,:,:],cmap=plt.cm.PuBu,levels=np.arange(100000,2050000,50000))
    cbar = plt.colorbar()
    cbar.set_label('total CO$_{2}O$ mixing ratio [ppmv]')
    plt.xlabel('distance [km]')
    plt.ylabel('height [m]')
    ax.set_yticks(np.arange(0,numLvl,5))
    ax.set_yticklabels(rx.lvl[::5])
    ax.set_xticks(np.arange(0,450,50))
    ax.set_xticklabels((np.arange(0,450,50)*0.040).astype(int))
    plt.title('EVOLUTION OF TOTAL CROSS-WIND CO$_2$')

    # update the frame
    def update_plot(n,cw_sum,cntr):
        del cntr
        cntr = plt.contourf(cw_sum[n,:,:],cmap=plt.cm.PuBu,levels=np.arange(100000,2050000,50000))
        return cntr,

    #plot all frames
    ani=animation.FuncAnimation(fig, update_plot, 300, fargs=(cw_sum,cntr), interval=1)
    ani.save(rx.fig_dir + 'Plume_CS_animation.mp4', writer='ffmpeg',fps=10)
    # plt.show()
