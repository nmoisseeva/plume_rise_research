##nmoisseeva@eoas.ubc.ca
#July 2019
#script to compare micromet-tower data  observed w wind speeds; plot total heat flux
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
from scipy import stats
import wrf
import imp
from datetime import datetime, timezone, timedelta
import matplotlib.dates as mdates
import pytz

#====================INPUT===================

#all common variables are stored separately
import rxcadreDRY as rx
imp.reload(rx) 	            #force load each time


datapath = './csv/W-H_micromet1Hz.csv'

freq = 1 		     #hz of wind measurements
ave_int_H = 300 	 #seconds
ave_int_W = 60       #seconds



hfx_sensor_conversion = 0.170           #mV m^2/kW
hfx_sensor_height = 2.8                 #meters AGL
w_sensor_height = 5.8                   #meters AGL
#======================end of input=======================

#get obs data
obs_data = np.genfromtxt(datapath,usecols=(1,2),skip_header=3,delimiter=',',dtype=float)
str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%y-%m-%d %H:%M') -  timedelta(hours=6)
time_data = np.genfromtxt(datapath,usecols=(0),skip_header=3,delimiter=',', converters = {0: str2date}, dtype=str)
obs_W = obs_data[:,0]
obs_H = obs_data[:,1]*1000/(hfx_sensor_conversion)

#get model data
print('Extracting NetCDF data from %s ' %rx.wrfdata)
nc_data = netcdf.netcdf_file(rx.wrfdata, mode ='r')
UTMx = nc_data.variables['XLONG'][0,:,:] + rx.ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + rx.ll_utm[1]

ncdict = wrf.extract_vars(nc_data, None, ('PHB','PH','W'))

#get height and destagger vars
print('...destaggering data')
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
w = wrf.destagger(ncdict['W'],1)

nT,nZ,nY,nX = np.shape(zstag)

#--------------NOT TESTED FROM HERE OWN-------------------

#create timeseries of micromet-tower wind
print('-->Creating OBS timeries')
num_pts = int(freq * ave_int_W)
data_samples = int(np.shape(obs_data)[0] / num_pts)
wMET = []
hMET = []
subset_time = time_data[::num_pts]
for nSample in range(data_samples):
    subset_W = obs_W[nSample*num_pts:(nSample+1)*num_pts]
    subset_H = obs_H[nSample*num_pts:(nSample+1)*num_pts]
    aveVal_W = np.mean(subset_W)
    aveVal_H = np.mean(subset_H)
    wMET.append(aveVal_W)
    hMET.append(aveVal_H)


#get grid location of where CSU is
print('-->Finding closest model location')
grid_coord_atm = np.array(list(zip(UTMy.ravel(),UTMx.ravel())))
gridTree = KDTree(grid_coord_atm)
METdist, METgrid_id = gridTree.query([rx.met_lcn[1], rx.met_lcn[0]])
METidx = np.unravel_index(METgrid_id,np.shape(UTMx))
#--------------NOT TESTED FROM HERE OWN-------------------

#get windspeed and height vector at the micromet tower location location
wrf_W = w[:,:,METidx[0],METidx[1]]
z_vector = np.mean(z[:,:,METidx[0],METidx[1]],0)

#create timeseries of WRF wind averaged to the same interval
wWRF = []
for nSample in range(rx.run_min):
	set_len = int(ave_int_W/rx.hist_int)
	subset = wrf_W[nSample*set_len:(nSample+1)*set_len,(0,1)]
	aveSet = np.mean(subset,0)
	wWRF.append(aveSet)


#=======================PLOTTING=======================
plt.title('MICROMET vs WRF  Vertical Velocity w ')
plt.plot(subset_time[-rx.run_min:],np.array(wWRF)[:,0], label='WRF-SFIRE 8m')
plt.plot(subset_time[-rx.run_min:],np.array(wMET)[-rx.run_min:], label='Micromet tower winds 5.8m')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('wind speed [m/s]')
plt.legend()
plt.show()
plt.close()

plt.title('MET-TOWER TOTAL HFX ')
plt.plot(subset_time[0:rx.run_min],np.array(hMET)[0:rx.run_min], label='Micromet tower Heat flux 2.8m')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('heat flux (W/m2)')
plt.legend()
plt.show()
plt.close()
