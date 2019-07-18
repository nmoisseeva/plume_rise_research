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
ave_int_W = 60       #min



hfx_sensor_conversion = 0.170           #mV m^2/kW
hfx_sensor_height = 2.8                 #meters AGL
w_sensor_height = 5.8                   #meters AGL
#======================end of input=======================

#get obs data
obs_data = np.genfromtxt(datapath,usecols=(1,2),skip_header=3,delimiter=',',dtype=float)
str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%y-%m-%d %H:%M') -  timedelta(hours=6)
time_data = np.genfromtxt(datapath,usecols=(0),skip_header=3,delimiter=',', converters = {0: str2date}, dtype=str)
obs_W = obs_data[:,0]
obs_H = obs_data[:,1]/hfx_sensor_conversion

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
num_pts = 60 * freq * ave_int
data_samples = np.shape(data)[0] / num_pts

#create timeseries of micromet-tower wind
print('-->Creating OBS timeries')
num_pts = int(freq * ave_int_W)
data_samples = int(np.shape(obs_data)[0] / num_pts)
uCSU = []
dirCSU = []
subset_time = time_data[::num_pts]
for nSample in range(data_samples):
	subset = obs_data[nSample*num_pts:(nSample+1)*num_pts,(0,2,4,6)]
	aveVal = np.mean(subset,0)
	uCSU.append(aveVal)

	dirset = obs_data[nSample*num_pts:(nSample+1)*num_pts,7]
	diraveVal = np.mean(dirset,0)
	dirCSU.append(diraveVal)
