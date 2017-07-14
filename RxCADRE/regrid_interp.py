# from scipy.io import netcdf
from Scientific.IO import NetCDF
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib as mpl 
from matplotlib import path 
from mpl_toolkits import basemap
import mpl_toolkits.basemap.pyproj as pyproj
import warnings
import os
warnings.filterwarnings("ignore")

datadir ='/Users/nmoisseeva/data/plume/RxCADRE'
savedir ='/regrid'
wrffile = 'wrfout_L2G_cat1'
# ll_utm = np.array([518800,3377000])
ll_utm = np.array([519500,3377000])

# ll_utm = np.array([523500,3377000])


#======================end of input=======================
print('REGRIDDING SCRIPT FOR WRF-IDEAL')

#copy netcdf file to a new directory
wrforig = datadir + '/' + wrffile
wrfnew = datadir + savedir + '/' + wrffile + '_regrid'
print('.....making a copy of NetCDF data in %s - will NOT overwrite!' %wrfnew)
cmnd = "cp -n %s %s" %(wrforig,wrfnew)
os.system(cmnd)

print('.....extracting NetCDF data from %s ' %wrfnew)
nc_data = NetCDF.NetCDFFile(wrfnew, 'a')

#create a UTM grid
UTMx = nc_data.variables['XLONG'][:,:,:] + ll_utm[0]
UTMy = nc_data.variables['XLAT'][:,:,:] + ll_utm[1]


#convert coordinate systems to something basemaps can read
print('.....transforming coordinates to WSG projection')
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())

WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))

nc_data.variables['XLONG'][:,:,:] = np.float32(WLONG)
nc_data.variables['XLAT'][:,:,:] = np.float32(WLAT)
nc_data.close()

