#nmoisseeva@eoas.ubc.ca
#Feb 2019
#script to compare CSU maps observed wind speeds with WRF windspeeds

import numpy as np
import matplotlib.pyplot as plt
from Scientific.IO import NetCDF
from scipy.io import netcdf
from scipy.spatial import KDTree
import wrf

datapath = './U.csv'
data = np.genfromtxt(datapath,usecols=(1,2,3,4,5,6,7,8),skip_header=3,delimiter=',',dtype=float)
csu_lcn = [525803.12, 3378544.15] #[30.539,86.731]
ll_utm = np.array([517000,3377000]) #Feb 2019

freq = 0.5 #hz of wind measurements
ave_int = 1 #min

wrfdata = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs'

#======================end of input=======================

print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')
UTMx = nc_data.variables['XLONG'][0,:,:] + ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + ll_utm[1]


ncdict = wrf.extract_vars(nc_data, None, ('PHB','PH','U','V'))

#get height and destagger vars
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
u = wrf.destagger(ncdict['U'],3)
v = wrf.destagger(ncdict['V'],2)

nT,nZ,nY,nX = np.shape(z)


#create timeseries of CSU wind
num_pts = int(60 * freq * ave_int)
data_samples = int(np.shape(data)[0] / num_pts)
uCSU = []
for nSample in range(data_samples):
	subset = data[nSample*num_pts:(nSample+1)*num_pts,(0,2,4,6)]
	ave1min = np.mean(subset,0)
	uCSU.append(ave1min)

#get grid location of where CSU is
grid_coord_atm = zip(UTMy.ravel(),UTMx.ravel())
gridTree = KDTree(grid_coord_atm)
CSUdist, CSUgrid_id = gridTree.query([csu_lcn[1], csu_lcn[0]])
CSUidx = np.unravel_index(CSUgrid_id,np.shape(UTMx))

#get windspeed and height vector at the CSU location
wrf_U = np.sqrt(u[:,:,CSUidx[0],CSUidx[1]]**2 +v[:,:,CSUidx[0],CSUidx[1]]**2)
z = np.mean(z[:,:,CSUidx[0],CSUidx[1]],0)

#!!!!HARDCODED!!!!!
#create timeseries of WRF wind averaged to the same interval (assumes 45min run and 10s history interval)
uWRF = []
for nSample in range(45):
	subset = wrf_U[nSample*6:(nSample+1)*6,(0,1)]
	print subset
	ave1min = np.mean(subset,0)
	uWRF.append(ave1min)

#=======================PLOTTING=======================
plt.title('CSU (6.2m,11.2m) vs WRF (~8m) WINDSPEEDS ')
plt.plot(np.array(uWRF)[:,0], label='WRF ~8m')
plt.plot(np.array(uCSU)[-45:,0], label='CSU winds 6.2m')
plt.plot(np.array(uCSU)[-45:,1], label='CSU winds 11.3m')
plt.xlabel('minutes')
plt.ylabel('wind speed [m/s]')
plt.legend()
plt.show()
plt.close()


plt.title('CSU (21.1m, 30.7m) vs WRF (~25m) WINDSPEEDS ')
plt.plot(np.array(uWRF)[:,1], label='WRF ~25m')
plt.plot(np.array(uCSU)[-45:,2], label='CSU winds 21.1m')
plt.plot(np.array(uCSU)[-45:,3], label='CSU winds 30.7m')
plt.xlabel('minutes')
plt.ylabel('wind speed [m/s]')
plt.legend()
plt.show()
plt.close()
