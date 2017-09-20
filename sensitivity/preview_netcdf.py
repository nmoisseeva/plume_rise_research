import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
from scipy.ndimage.interpolation import zoom
import matplotlib.animation as animation
from matplotlib import path 
from mpl_toolkits import basemap
import mpl_toolkits.basemap.pyproj as pyproj
import os.path
import pickle
import mpl_toolkits.mplot3d.axes3d as p3
import mpl_toolkits.mplot3d as a3
from matplotlib import animation

#====================INPUT===================
wrfdata = '/Users/nmoisseeva/data/plume/sensitivity/raw/wrfout_W6S400F3R0'

#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  

tstep = round((nc_data.variables['XTIME'][1]- nc_data.variables['XTIME'][0]) * 60.)  #timestep in sec
# WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]

# ghfx = np.copy(nc_data.variables['GRNHFX'][:,:,:]) 	#extract fire heat flux
T = np.copy(nc_data.variables['T'][80,:,:,:])
Tprofile = T[:,100,100]
plt.plot(Tprofile)
plt.plot(nc_data.variables['T'][10,:,100,100])
plt.plot(nc_data.variables['T'][50,:,100,100])
plt.plot(nc_data.variables['T'][0,:,100,100])
plt.show()




Tprofile = np.nanmean(np.nanmean(T,2),1)



z = (nc_data.variables['PHB'][0,:,100,100] + nc_data.variables['PH'][0,:,100,100])/9.81