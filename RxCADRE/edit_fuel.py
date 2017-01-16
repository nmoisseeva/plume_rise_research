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

# wrfinput='/Users/nadya2/Applications/WRFV3/test/em_fire/wrfinput_d01'
wrfinput='/Users/nadya2/Applications/WRF-SFIRE/wrf-fire/WRFV3/test/em_fire/rxcadre/wrfinput_d01'


# ll_utm = np.array([521620,3376766]) - first domain
# ll_utm = np.array([518317,3373798])
# ll_utm = np.array([520317,3375798])
ll_utm = np.array([518800,3377000])


bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'

#four lines (strip headfire method) walking ignition
fire_dict_utm = {'fireline1':{'start':np.array([525828,3379011]), 'end':np.array([524551,3378179])},\
				'fireline2':{'start':np.array([525729,3379075]), 'end':np.array([524487,3378275])},\
				'fireline3':{'start':np.array([525612,3379181]), 'end':np.array([524409,3378388])},\
				'fireline4':{'start':np.array([525538,3379244]), 'end':np.array([524331,3378480])} }
fuel_cat = 3

#======================end of input=======================
print('Extracting NetCDF data from %s ' %wrfinput)
nc_data = NetCDF.NetCDFFile(wrfinput, 'a')

#create a UTM grid
UTMx = nc_data.variables['XLONG'][0,:,:] + ll_utm[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + ll_utm[1]
UTMfx = nc_data.variables['FXLONG'][0,:,:] + ll_utm[0]
UTMfy = nc_data.variables['FXLAT'][0,:,:] + ll_utm[1]

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
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)

print('..... replacing fuel data')
l2g = path.Path(bm.fire_bounds[5])
l2g_mask = l2g.contains_points(zip(WGSfx,WGSfy))
l2g_mask = np.reshape(l2g_mask, np.shape(UTMfx))
fuel = nc_data.variables['NFUEL_CAT'][0,:,:]
fuel[l2g_mask] = fuel_cat
fuel[~l2g_mask] = 14

bm.contourf(WLONGf, WLATf,fuel)
plt.colorbar()
plt.show()

nc_data.variables['NFUEL_CAT'][0,:,:] = fuel
nc_data.close()

#convert fireline cooredinate to WRF input
print('..... converting fireline cooredinates to WRF input')
for key in fire_dict_utm:
	fstart = fire_dict_utm[key]['start'][:] - ll_utm[:]
	fend = fire_dict_utm[key]['end'][:] - ll_utm[:]
	print'     ..... %s: start (%sm, %sm) - end (%sm, %sm)' %(key, fstart[0],fstart[1],fend[0],fend[1])

