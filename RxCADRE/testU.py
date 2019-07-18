#nmoisseeva@eoas.ubc.ca
#Feb 2019
#script to compare CSU maps observed wind speeds with WRF windspeeds

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
import rxcadreMOIST as rx
imp.reload(rx) 	            #force load each time


datapath = './csv/WSPD-CSU-MAPS.csv'

freq = 0.5 		#hz of wind measurements
ave_int = 60 	#seconds
#======================end of input=======================

#get obs data
obs_data = np.genfromtxt(datapath,usecols=(1,2,3,4,5,6,7,8),skip_header=3,delimiter=',',dtype=float)
str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%y-%m-%d %H:%M') -  timedelta(hours=6)
time_data = np.genfromtxt(datapath,usecols=(0),skip_header=3,delimiter=',', converters = {0: str2date}, dtype=str)

#get model data
print('Extracting NetCDF data from %s ' %rx.wrfdata)
nc_data = netcdf.netcdf_file(rx.wrfdata, mode ='r')
UTMx = nc_data.variables['XLONG'][0,:,:] + rx.ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + rx.ll_utm[1]

ncdict = wrf.extract_vars(nc_data, None, ('PHB','PH','U','V'))

#get height and destagger vars
print('...destaggering data')
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)

nT,nZ,nY,nX = np.shape(z)

#create timeseries of CSU wind
print('-->Creating OBS timeries')
num_pts = int(freq * ave_int)
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

#get grid location of where CSU is
print('-->Finding closest model location')
grid_coord_atm = np.array(list(zip(UTMy.ravel(),UTMx.ravel())))
gridTree = KDTree(grid_coord_atm)
CSUdist, CSUgrid_id = gridTree.query([rx.csu_lcn[1], rx.csu_lcn[0]])
CSUidx = np.unravel_index(CSUgrid_id,np.shape(UTMx))

#get windspeed and height vector at the CSU location
wrf_U = np.sqrt(u[:,:,CSUidx[0],CSUidx[1]]**2 +v[:,:,CSUidx[0],CSUidx[1]]**2)
z = np.mean(z[:,:,CSUidx[0],CSUidx[1]],0)

#------------WHY IS 1MIN HARDCODED????-----------------
#create timeseries of WRF wind averaged to the same interval (assumes 45min run and 10s history interval) -
uWRF = []
for nSample in range(rx.run_min):
	set_len = int(60/rx.hist_int)
	subset = wrf_U[nSample*set_len:(nSample+1)*set_len,(0,1)]
	ave1min = np.mean(subset,0)
	uWRF.append(ave1min)

#=======================PLOTTING=======================
plt.title('CSU (6.2m,11.2m) vs WRF (8m) WINDSPEEDS ')
plt.plot(subset_time[-rx.run_min:],np.array(uWRF)[:,0], label='WRF-SFIRE 8m')
plt.plot(subset_time[-rx.run_min:],np.array(uCSU)[-rx.run_min:,0], label='CSU winds 6.2m')
plt.plot(subset_time[-rx.run_min:],np.array(uCSU)[-rx.run_min:,1], label='CSU winds 11.3m')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('wind speed [m/s]')
plt.legend()
plt.show()
plt.close()


plt.title('CSU (21.1m, 30.7m) vs WRF (25m) WINDSPEEDS ')
plt.plot(subset_time[-rx.run_min:],np.array(uWRF)[:,1], label='WRF-SFIRE 25m')
plt.plot(subset_time[-rx.run_min:],np.array(uCSU)[-rx.run_min:,2], label='CSU winds 21.1m')
plt.plot(subset_time[-rx.run_min:],np.array(uCSU)[-rx.run_min:,3], label='CSU winds 30.7m')
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('wind speed [m/s]')
plt.legend()
plt.show()
plt.close()



#------create a regression line for wind------

xaxis = np.arange(122)
m, b, r_value, p_value, std_err = stats.linregress(xaxis,np.array(dirCSU)[-122:])
fit = b + m * xaxis

plt.title('CSU (30.7m) WIND DIRECTION')
plt.scatter(xaxis,np.array(dirCSU)[-122:])
plt.plot(xaxis,fit,'r--')
plt.xlabel('time (CST)')
plt.ylabel('wind direction [deg]')
plt.ylim([40,220])
ax = plt.gca()
ax.set_xticks(np.arange(0,121,30))
ax.set_xticklabels(['11:10','11:40','12:10','12:40','13:10'])
plt.savefig(rx.fig_dir + 'WindDir.pdf')
plt.show()
plt.close()
