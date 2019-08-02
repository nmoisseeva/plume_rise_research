#April 2019
#nmoisseeva@eoas.ubc.ca

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import wrf
import datetime as dt
import matplotlib.dates as mdates
import pandas as pd
from matplotlib.ticker import MaxNLocator
import imp
import sys
from scipy.spatial import KDTree


#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx)		        #force load each time

#-----------------edited from within steps below------

basetime = dt.datetime(year=2012,month=11,day=10)

# #load model data
# qinterp = np.load(rx.interp_path)   # load here the above pickle
# print('Interpolated data found at: %s' %rx.interp_path)

#extract and format dispersion data
print('Importing dispersion data from %s' %rx.disp_data)
disp_dict = {}
disp_array = np.genfromtxt(rx.disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')
disp_dict['time'] = np.array([basetime + dt.timedelta(seconds = int(i)) for i in disp_array[:,0]])
disp_dict['CO'] = disp_array[:,1]
disp_dict['CO2'] = disp_array[:,2]
disp_dict['CH4'] = disp_array[:,3]
disp_dict['H2O'] = disp_array[:,4]
disp_dict['lcn'] = np.array(list(zip(disp_array[:,5],disp_array[:,6],disp_array[:,7]-rx.sfc_hgt)))
disp_dict['meta']= 'time: timestamp| \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'

prof_dt_start = [dt.datetime.combine(basetime,dt.datetime.strptime(i, '%H:%M:%S').time()) for i in rx.profiles_start]
prof_dt_end = [dt.datetime.combine(basetime,dt.datetime.strptime(i, '%H:%M:%S').time()) for i in rx.profiles_end]
cs_start, cs_end = dt.datetime.combine(basetime,dt.datetime.strptime(rx.corkscrew[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(rx.corkscrew[1], '%H:%M:%S').time())
gg_start, gg_end = dt.datetime.combine(basetime,dt.datetime.strptime(rx.garage[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(rx.garage[1], '%H:%M:%S').time())

#import wrf data
print('Extracting NetCDF data from %s ' %rx.spinup_path)
wrfdata = netcdf.netcdf_file(rx.spinup_path, mode ='r')
ncdict = wrf.extract_vars(wrfdata, None, ('T','PHB','PH','QVAPOR','XTIME'))

#get geopotential array and convert to height
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
z_ave = np.mean(z, (0,2,3))
# # START AT PATH EXTRACTION (FOR HIGH FREQUENCY DATA)
#
# model_datetime = [basetime + dt.timedelta(hours=int(rx.runstart[:2]),minutes=t) for t in ncdict['XTIME']]
# #get indecies of samples corresponding to model output times
# tidx = [np.argmin(abs(disp_dict['time'] - t)) for t in model_datetime] #times since start
#
# #construct KDtree from idealized grid
# lat = wrfdata.variables['XLAT'][0,:,:]
# lon = wrfdata.variables['XLONG'][0,:,:]
# grid_coord = list(zip(lat.ravel(),lon.ravel()))
# gridTree = KDTree(grid_coord)
# dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx][:,0:2])
#
#
# # mod_val, obs_val = np.empty((len(tidx))) * np.nan, np.empty((len(tidx))) * np.nan
# # obs_h = []
# # obs_lon, obs_lat = [],[]
# # print('Finding nearest points....')
# # for nt in range(len(tsec)):
# #     print('...tstep: %s') %nt
# #     idxy,idxx = np.unravel_index(grid_id[nt],np.shape(lat))
# #     idxz = np.argmin(abs(lvl-disp_dict['lcn'][tidx][nt][2]))
# #     mod_val[nt] = qinterp[nt,idxz,idxy,idxx]
# #     obs_val[nt] = disp_dict['CO2'][tidx[nt]]		#in percent by volume
# #     obs_h.append(disp_dict['lcn'][tidx][nt][2])
# #     obs_lon.append(disp_dict['lcn'][tidx][nt][1])
# #     obs_lat.append(disp_dict['lcn'][tidx][nt][0])


#===============================PLOTTING===========================
#plot plame height vs time
plt.figure()
plt.title('AIRPLANE HEIGHT')
plt.plot(disp_dict['time'],disp_dict['lcn'][:,2])
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
for nProf in range(len(prof_dt_end)):
	ax.fill_between(pd.date_range(prof_dt_start[nProf],prof_dt_end[nProf],freq='S'),0,2000,facecolor='green', alpha=0.2)
ax.fill_between(pd.date_range(cs_start,cs_end,freq='S'),0,2000,facecolor='red', alpha=0.2)
ax.fill_between(pd.date_range(gg_start,gg_end,freq='S'),0,2000,facecolor='red', alpha=0.2)
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('airplane height [m]')
plt.ylim([0,2000])
plt.savefig(rx.fig_dir + 'PlaneHeightProfiles.pdf')
plt.show()


#CONSIDER EXtRACTING AVERAGE LOCATION OF AIRPLAIN PROFILE AND GETTING THE DATA FROM THERE
if rx.moist_run:
    fig = plt.figure(figsize=(16,5))
    lbl = ['(a)','(b)','(c)']
    for nProfile in range(2):
        model_min = (prof_dt_end[nProfile].hour - int(rx.runstart[:2]))*60 + prof_dt_end[nProfile].minute
        model_time_idx = np.argmin(abs(ncdict['XTIME'] - model_min))
        print("Model timing off from observation by: %s minutes" %min(abs(ncdict['XTIME'] - model_min)))
        i = np.argmin(abs(disp_dict['time'] - prof_dt_start[nProfile]))
        f = np.argmin(abs(disp_dict['time'] - prof_dt_end[nProfile]))
        ax = plt.subplot(1,3,nProfile+1)
        plt.ylabel('hight AGL [m]')
        ax = plt.gca()
        model_mean = np.mean(ncdict['QVAPOR'][model_time_idx,:,:,:],(1,2))
        print(model_mean)
        model_vmr =  (28.9644 / 18.01528) * model_mean           #convert to percent by volume molar mass(dry air/h20)
        plt.plot(model_vmr*100,z_ave,color='tab:red',label='model water vapour'  )
        plt.plot(disp_dict['H2O'][i:f],disp_dict['lcn'][i:f,2],color='tab:blue',label='$H_20$ profile from flight')
        plt.legend()
        # ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel('$H_2O$ mixing ratio [%]')
    plt.tight_layout()
    plt.savefig(rx.fig_dir + 'BLGrowthEvaluation.pdf')
    plt.show()
else:
    fig = plt.figure(figsize=(16,5))
    lbl = ['(a)','(b)','(c)']
    for nProfile in range(3):
        model_min = (prof_dt_end[nProfile].hour - int(rx.runstart[:2]))*60 + prof_dt_end[nProfile].minute
        i = np.argmin(abs(disp_dict['time'] - prof_dt_start[nProfile]))
        f = np.argmin(abs(disp_dict['time'] - prof_dt_end[nProfile]))
        ax = plt.subplot(1,3,nProfile+1)
        plt.ylabel('hight AGL [m]')
        ax1 = plt.gca()
        color = 'tab:red'
        ax1.set_xlabel('height [m]')
        ax1.set_xlabel('potential temperature [K]', color=color)
        l1 = ax1.plot(ncdict['T'][model_min,:,0,0]+300,z_ave,color=color,label='model potential temperature' )
        ax1.set_xlim([290,310])
        ax1.tick_params(axis='x', labelcolor=color)
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.legend()
        ax2 = ax1.twiny()
        color = 'tab:blue'
        ax2.set_xlabel('$H_2O$ mixing ratio [%]', color=color)
        l2 = ax2.plot(disp_dict['H2O'][i:f],disp_dict['lcn'][i:f,2],color=color,label='$H_20$ profile from flight')
        ax2.set_xlim([1.2,0.2])
        ax2.tick_params(axis='x', labelcolor=color)
        lns = l1 + l2
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc=2,title='%s %s CST' %(lbl[nProfile],rx.profiles_start[nProfile]))
    plt.tight_layout()
    plt.savefig(rx.fig_dir + 'BLGrowthEvaluation.pdf')
    plt.show()
