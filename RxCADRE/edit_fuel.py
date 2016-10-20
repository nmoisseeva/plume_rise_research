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

wrfinput='/Users/nadya2/Applications/WRFV3/test/em_fire/wrfinput_d01'

ll_utm = np.array([521620,3376766])
bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'

fire_dict_utm = {'fireline1':{'start':np.array([525823,3379005]), 'end':np.array([525151,3378271])},\
				'fireline2':{'start':np.array([525625,3379182]), 'end':np.array([525005,3378787])},\
				'fireline3':{'start':np.array([525151,3378571]), 'end':np.array([524564,3378184])},\
				'fireline4':{'start':np.array([525005,3378787]), 'end':np.array([524415,3378389])} }

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
# import shape file
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)
plt.show()


print('..... replacing fuel data')
l2g = path.Path(bm.fire_bounds[5])
l2g_mask = l2g.contains_points(zip(WGSfx,WGSfy))
l2g_mask = np.reshape(l2g_mask, np.shape(UTMfx))
fuel = nc_data.variables['NFUEL_CAT'][0,:,:]
fuel[l2g_mask] = 3
fuel[fuel!=3] = 14

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



# plt.figure(figsize=(6,6))
# plt.plot(x[:,1], x[:,0])
# plt.xlim([295,315])
# plt.ylim([0,2000])
# plt.ylabel('height [m]')
# plt.xlabel('potential temperature [K]')
# plt.tight_layout()
# plt.savefig('/Users/nadya2/GoogleDrive/PhD/comps_written/roland/figs/sounding.pdf')