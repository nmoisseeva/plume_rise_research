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

#====================INPUT===================
wrfpath = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs_spinup'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile
sfc_hgt = 62 							#surface height MSL (m)
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)
#-------------------end of input--------------

# profiles_sec = np.array([[7367,7507],[7985,8357],[8813,8982]])
profiles_start = ['12:02:30','12:13:00','12:27:00']
profiles_end = ['12:05:30','12:19:00','12:29:30']
garage = ['12:36:00','13:04:00']
corkscrew = ['13:06:30','13:12:00']

#-----------------edited from within steps below------

model_ssm = int(runstart[0:2])*3600 + int(runstart[3:5])*60

#load model data
qinterp = np.load(interp_path)   # load here the above pickle
print('Interpolated data found at: %s' %interp_path)

#extract and format dispersion data
print('Importing dispersion data from %s' %disp_data)
disp_dict = {}
disp_array = np.genfromtxt(disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')
start_idx = np.argmin(abs(disp_array[:,0] - model_ssm))  #find simulation start time index
disp_dict['time']= (disp_array[start_idx:,0] - model_ssm +1)
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]-sfc_hgt))
disp_dict['meta']= 'time: min since restart run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'

#convert all timestamps to datetime objects
basetime = dt.datetime(year=2012,month=11,day=10)
timestamp = [basetime + dt.timedelta(hours = 10, seconds = i) for i in disp_dict['time']]

prof_dt_start = [dt.datetime.combine(basetime,dt.datetime.strptime(i, '%H:%M:%S').time()) for i in profiles_start]
prof_dt_end = [dt.datetime.combine(basetime,dt.datetime.strptime(i, '%H:%M:%S').time()) for i in profiles_end]
cs_start, cs_end = dt.datetime.combine(basetime,dt.datetime.strptime(corkscrew[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(corkscrew[1], '%H:%M:%S').time())
gg_start, gg_end = dt.datetime.combine(basetime,dt.datetime.strptime(garage[0], '%H:%M:%S').time()),dt.datetime.combine(basetime,dt.datetime.strptime(garage[1], '%H:%M:%S').time())

#load pre-burn moisture profile
pre_burn_qv = np.genfromtxt(pre_moisture, skip_header=1, delimiter=',')

#import wrf data
print('Extracting NetCDF data from %s ' %wrfpath)
wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')
ncdict = wrf.extract_vars(wrfdata, None, ('T','PHB','PH'))

#get geopotential array and convert to height
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
z_ave = np.mean(z, (0,2,3))

#===============================PLOTTING===========================
#plot plame height vs time
plt.figure()
plt.title('AIRPLANE HEIGHT')
plt.plot(timestamp,disp_dict['lcn'][:,2])
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
plt.savefig(fig_dir + 'PlaneHeightProfiles.pdf')
plt.show()



# #create subplots of profiles
# fig = plt.figure(figsize=(16,12))
# # ax = plt.subplot(2,2,1)
# color = 'tab:red'
# ax1 = plt.gca()
# ax1.set_xlabel('height [m]')
# ax1.set_xlabel('potential temperature [K]', color=color)
# l1 = ax1.plot(ncdict['T'][0,:,0,0]+300,z_ave,color=color,label='intial model $T$ from 10:00:00 CST sounding')
# ax1.set_xlim([291.5,310])
# ax1.tick_params(axis='x', labelcolor=color)
# ax2 = ax1.twiny()
# color = 'tab:blue'
# ax2.set_xlabel('$H_2O$ mixing ratio [%]', color=color)
# l2 = ax2.plot(pre_burn_qv[:,1],pre_burn_qv[:,0],color=color,label='$H_20$ profile from 10:00:00 CST sounding')
# ax2.set_xlim([1.2,0.2])
# ax2.tick_params(axis='x', labelcolor=color)
# lns = l1 + l2 								#combine lines to create a single legend
# labs = [l.get_label() for l in lns] 		#get labels
# ax.legend(lns, labs, loc=2)

fig = plt.figure(figsize=(16,5))
lbl = ['(a)','(b)','(c)']
for nProfile in range(3):
    model_time = profiles_start[nProfile]
    i, f = profiles_sec[nProfile,:]
    i = np.argmin(abs(disp_dict['time'] - i))
    f = np.argmin(abs(disp_dict['time'] - f))
    ax = plt.subplot(1,3,nProfile+1)
    plt.ylabel('hight AGL [m]')
    ax1 = plt.gca()
    color = 'tab:red'
    ax1.set_xlabel('height [m]')
    ax1.set_xlabel('potential temperature [K]', color=color)
    l1 = ax1.plot(ncdict['T'][model_min,:,0,0]+300,z_ave,color=color,label='model $T$' )
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
    ax.legend(lns, labs, loc=2,title='%s %s CST' %(lbl[nProfile],model_time))
plt.tight_layout()
plt.savefig(fig_dir + 'BLGrowthEvaluation.pdf')
plt.show()
