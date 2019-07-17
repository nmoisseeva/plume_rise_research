from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib as mpl
from matplotlib import path
from mpl_toolkits import basemap
import pyproj as pyproj
import warnings
import os
warnings.filterwarnings("ignore")
import imp

#====================INPUT===================

#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx) 	            #force load each time

wrfinput='/Users/nmoisseeva/sfire/wrf-fire/WRFV3/test/em_fire/rxcadre_moist/wrfinput_d01'
input_fc = '/Users/nmoisseeva/sfire/wrf-fire/WRFV3/test/em_fire/rxcadre_moist/input_fc'

#======================end of input=======================
print('Extracting NetCDF data from %s ' %wrfinput)
nc_data = netcdf.netcdf_file(wrfinput, mode ='a')

#create a UTM grid
UTMx = nc_data.variables['XLONG'][0,:,:] + rx.ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + rx.ll_utm[1]
UTMfx = nc_data.variables['FXLONG'][0,:,:] + rx.ll_utm[0]
UTMfy = nc_data.variables['FXLAT'][0,:,:] + rx.ll_utm[1]

#convert coordinate systems to something basemaps can read
print('..... transforming coordinates to WSG projection')
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())
WGSfx, WGSfy= pyproj.transform(epsg26916,wgs84,UTMfx.ravel(),UTMfy.ravel())

WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))
WLONGf, WLATf = np.reshape(WGSfx, np.shape(UTMfx)), np.reshape(WGSfy, np.shape(UTMfy))

#generate basemap
print('..... configuring basemaps')
bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)

#pull landsat image and save for future use (fabour)
landsat_pull_cmnd = "bash getWMSImage.sh -o ./npy/landsat_%s_%s_%s_%s.tiff -m landsat %s %s %s %s" \
				%(WLONG[0,0],WLAT[0,0],WLONG[-1,-1],WLAT[-1,-1],WLONG[0,0],WLAT[0,0],WLONG[-1,-1],WLAT[-1,-1])
os.system(landsat_pull_cmnd)

#add plot shapefile
polygons = bm.readshapefile(rx.bounds_shape,name='fire_bounds',drawbounds=True)

print('..... replacing fuel data')
l2g = path.Path(bm.fire_bounds[5])
print('.......-select points')
l2g_mask = l2g.contains_points(zip(WGSfx,WGSfy))
l2g_mask = np.reshape(l2g_mask, np.shape(UTMfx))

print('.......-copying fuel data')
fuel = nc_data.variables['NFUEL_CAT'][0,:,:]
fuel[l2g_mask] = rx.fuel_cat
fuel[~l2g_mask] = 14

np.savetxt(input_fc, fuel.T,delimiter=' ',fmt='%d', header = '%s,%s' %(nc_data.dimensions['west_east_subgrid'],nc_data.dimensions['south_north_subgrid']))

#sanity-check plot
bm.contourf(WLONGf, WLATf,fuel)
plt.colorbar()
plt.show()

nc_data.close()

#convert fireline cooredinate to WRF input
print('..... converting fireline cooredinates to WRF input')
for key in rx.fire_dict_utm:
	fstart = rx.fire_dict_utm[key]['start'][:] - rx.ll_utm[:]
	fend = rx.fire_dict_utm[key]['end'][:] - rx.ll_utm[:]
	print('..... %s: start (%sm, %sm) - end (%sm, %sm)' %(key, fstart[0],fstart[1],fend[0],fend[1]))
