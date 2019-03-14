
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import wrf



#====================INPUT===================
wrfpath = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs_spinup'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile
sfc_hgt = 62 							#surface height MSL (m)
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)
#-------------------end of input--------------

profiles_sec = np.array([[7367,7507],[7985,8357],[8813,8982]])


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
disp_dict['time']= disp_array[start_idx:,0] - model_ssm +1
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]-sfc_hgt))
disp_dict['meta']= 'time: seconds since restart run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'

#plot plame height vs time
plt.figure()
plt.plot(disp_dict['time'],disp_dict['lcn'][:,2])
plt.show()

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


plt.figure(figsize=(16,12))
plt.subplot(2,2,1)
color = 'tab:red'
ax1 = plt.gca()
ax1.set_xlabel('height [m]')
ax1.set_xlabel('perturbation temperature [K]', color=color)
ax1.plot(ncdict['T'][0,:,0,0],z_ave,color=color,label='model profile 10am')
ax1.set_xlim([-8.5,10])
ax1.tick_params(axis='x', labelcolor=color)
ax2 = ax1.twiny()
color = 'tab:blue'
ax2.set_xlabel('water vapour (x-1)', color=color)
ax2.plot(pre_burn_qv[:,1]*-1,pre_burn_qv[:,0],color=color,label='sounding data 10am')
ax2.set_xlim([-1.2,-0.2])
ax2.tick_params(axis='x', labelcolor=color)


for nProfile in range(3):
    model_min = int(profiles_sec[nProfile,0]/60)
    print model_min
    i, f = profiles_sec[nProfile,:]
    i = np.argmin(abs(disp_dict['time'] - i))
    f = np.argmin(abs(disp_dict['time'] - f))
    plt.subplot(2,2,nProfile+2)
    ax1 = plt.gca()
    color = 'tab:red'
    ax1.set_xlabel('height [m]')
    ax1.set_xlabel('perturbation temperature [K]', color=color)
    ax1.plot(ncdict['T'][model_min,:,0,0],z_ave,color=color,label='model BL %s min' %model_min)
    ax1.set_xlim([-8.5,10])
    ax1.tick_params(axis='x', labelcolor=color)
    plt.legend()
    ax2 = ax1.twiny()
    color = 'tab:blue'
    ax2.set_xlabel('water vapour (*-1)', color=color)
    ax2.plot(disp_dict['H2O'][i:f]*-1,disp_dict['lcn'][i:f,2],color=color,label='flight data %s min' %model_min)
    ax2.set_xlim([-1.2,-0.2])
    ax2.tick_params(axis='x', labelcolor=color)
    plt.legend()
plt.show()
