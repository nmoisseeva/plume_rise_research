#This code compares wrf-spinup profiles with rxcadre radiometer profiler data
#jul 2019
#nmoisseeva@eoas.ubc.ca

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from datetime import datetime, timedelta
import wrf
import imp


#====================INPUT===================
#all common variables are stored separately
import rxcadreDRY as rx
imp.reload(rx) 	            #force load each time

#-----------------end of input------


print('===============RADIOMETER PROFILER DATA COMPARISON============')
print('Extracting observational data from %s ' %rx.radiometer_data)
obs_data = np.genfromtxt(rx.radiometer_data,skip_header=3,delimiter=',',dtype=float)
str2date = lambda x: datetime.strptime(x.decode('utf-8'), '%y/%m/%d %H:%M') -  timedelta(hours=6)
time_data = np.genfromtxt(rx.radiometer_data,usecols=(0),skip_header=3,delimiter=',', converters = {0: str2date}, dtype=str)[1:]

#clean up and sort obs data
profile_data = obs_data[:,2:]           #remove time and title rows
height_levels = obs_data[0,:] * 1000    #first row is height level in km, convert to meters

#import model spinup data
print('Extracting NetCDF data from %s ' %rx.spinup_path)
nc_data = netcdf.netcdf_file(rx.spinup_path, mode ='r')
UTMx = nc_data.variables['XLONG'][0,:,:] + rx.ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + rx.ll_utm[1]

ncdict = wrf.extract_vars(nc_data, None, ('PHB','PH','T'))

#get height and destagger vars
print('...destaggering data')
zstag = (ncdict['PHB'] + ncdict['PH'])/ 9.81
z = wrf.destagger(zstag,1)
t = ncdict['T'][:] + 300
nT,nZ,nY,nX = np.shape(z)

#get grid location of where the radiometer is
print('...inding closest model location')
grid_coord_atm = np.array(list(zip(UTMy.ravel(),UTMx.ravel())))
gridTree = KDTree(grid_coord_atm)
RADdist, RADgrid_id = gridTree.query([rx.rad_lcn[1], rx.rad_lcn[0]])
RADidx = np.unravel_index(RADgrid_id,np.shape(UTMx))

#get temperature profiles at radiometer location
wrf_T = t[:,:,RADidx[0],RADidx[1]]
z = np.mean(z[:,:,RADidx[0],RADidx[1]],0)


#---------------NEEED TO TEST CODE AND MATCH TIMING NEXT-----
