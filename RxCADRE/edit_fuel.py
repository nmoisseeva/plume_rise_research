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
warnings.filterwarnings("ignore")

# wrfinput='/Users/nadya2/Applications/WRFV3/test/em_fire/wrfinput_d01'
wrfinput='/Users/nadya2/Applications/WRF-SFIRE/wrf-fire/WRFV3/test/em_fire/rxcadre/wrfinput_d01'


# ll_utm = np.array([521620,3376766]) - first domain
# ll_utm = np.array([518317,3373798])
ll_utm = np.array([520317,3375798])

bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'

# #fourlines 1/2 length of plot (every other real ignition line, changed from 1/3 length each time)
# fire_dict_utm = {'fireline1':{'start':np.array([525813,3379004]), 'end':np.array([525105,3378544])},\
# 				'fireline2':{'start':np.array([525614,3379181]), 'end':np.array([524962,3378764])},\
# 				'fireline3':{'start':np.array([525115,3378546]), 'end':np.array([524565,3378184])},\
# 				'fireline4':{'start':np.array([524925,3378725]), 'end':np.array([524421,3378389])} }

#four lines (strip headfire method) walking ignition
fire_dict_utm = {'fireline1':{'start':np.array([525828,3379011]), 'end':np.array([524551,3378179])},\
				'fireline2':{'start':np.array([525729,3379075]), 'end':np.array([524487,3378275])},\
				'fireline3':{'start':np.array([525612,3379181]), 'end':np.array([524409,3378388])},\
				'fireline4':{'start':np.array([525538,3379244]), 'end':np.array([524331,3378480])} }
fuel_cat = 3

#======================end of input=======================
print('Extracting NetCDF data from %s ' %wrfinput)
nc_data = NetCDF.NetCDFFile(wrfinput, 'a')


#get dimensions of the data
nT,nY,nX = np.shape(nc_data.variables['XLONG'])

#construct grid
grid_coord = zip(nc_data.variables['XLONG'][0].ravel(),nc_data.variables['XLAT'][0].ravel())
fire_grid_coord = zip(nc_data.variables['FXLONG'][0].ravel(),nc_data.variables['FXLAT'][0].ravel())

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

#generate basemap
print('..... configuring basemaps')
bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)

#add plot shapefile
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)

print('..... replacing fuel data')
l2g = path.Path(bm.fire_bounds[5])
l2g_mask = l2g.contains_points(zip(WGSfx,WGSfy))
l2g_mask = np.reshape(l2g_mask, np.shape(UTMfx))
fuel = nc_data.variables['NFUEL_CAT'][0,:,:]
fuel[l2g_mask] = fuel_cat
fuel[~l2g_mask] = 14
# fuel[fuel!=fuel_cat] = 14


plt.contourf(fuel)
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

