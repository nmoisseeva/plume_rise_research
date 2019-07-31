# from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
# import matplotlib as mpl
# from matplotlib import path
from mpl_toolkits import basemap
import pyproj
import imp
import pickle
import os.path
#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx)		        #force load each time

#======================end of input=======================
print('GEOREFERENCING SCRIPT FOR WRF-IDEAL')

print('.....extracting NetCDF data from %s ' %rx.spinup_path)
ncdata = netcdf.netcdf_file(rx.spinup_path, mode ='r')

#create a UTM grid
UTMx = ncdata.variables['XLONG'][0,:,:] + rx.ll_utm[0]
UTMy = ncdata.variables['XLAT'][0,:,:] + rx.ll_utm[1]
FUTMx = ncdata.variables['FXLONG'][0,:,:] + rx.ll_utm[0]
FUTMy = ncdata.variables['FXLAT'][0,:,:] + rx.ll_utm[1]

#convert coordinate systems to something basemaps can read
print('.....transforming coordinates to WSG projection')
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())
FWGSx, FWGSy= pyproj.transform(epsg26916,wgs84,FUTMx.ravel(),FUTMy.ravel())

WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))
FWLONG, FWLAT = np.reshape(FWGSx, np.shape(FUTMx)), np.reshape(FWGSy, np.shape(FUTMy))

#save data in a dictionary
wrfgeo = {'XLONG':UTMx, 'XLAT':UTMy, 'FXLONG':FUTMx, 'FXLAT':FUTMy, 'WLONG':WLONG, 'WLAT':WLAT, 'FWLONG':FWLONG, 'FWLAT':FWLAT, \
            'metadata': 'XLONG, XLAT, FXLONG, FXLAT: amospheric and fire grid in UTM coordinates EPSG26916 (i.e. model output in m + ll._utm) | \
                        WLONG, WLAT, FWLONG, FWLAT: atmospheric and fire grid in WG84 coordinates for basemaps (EPSG4326)'}

np.save(rx.geo_path, wrfgeo)


#open/generate basemap
if os.path.isfile(rx.basemap_path):
	bm = pickle.load(open(rx.basemap_path,'rb'))   # load here the above pickle
	print('.....domain basemap found at: %s' %rx.basemap_path)
else:
	print('WARNING: no existing basemaps found: configuring a new basemap')
	bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)
	pickle.dump(bm,open(rx.basemap_path,'wb'),-1)  	# pickle the new map for later
	print('.....new basemap instance saved as: %s' %rx.basemap_path)


ncdata.close()
