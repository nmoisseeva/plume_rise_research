# referencing simulation with flight data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
import matplotlib.animation as animation
from matplotlib import path 
from mpl_toolkits import basemap
import mpl_toolkits.basemap.pyproj as pyproj


#====================INPUT===================
fireline1_utm = (525812,3378990)
fireline1_xy = (3400,1000)
wrfdata = '/Users/nadya2/data/plume/RxCADRE/wrfout'
bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'




#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  

#get dimensions of the data
nT,nY,nX = np.shape(nc_data.variables['XLONG'])

#construct KDtree from idealized grid
grid_coord = zip(nc_data.variables['XLONG'][0].ravel(),nc_data.variables['XLAT'][0].ravel())
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(fireline1_xy)

#create a UTM grid
UTMx = nc_data.variables['XLONG'][0,:,:] + fireline1_utm[0] - fireline1_xy[0]
UTMy = nc_data.variables['XLAT'][0,:,:] + fireline1_utm[1] - fireline1_xy[1]

#convert coordinate systems to something basemaps can read
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())

WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))

#generate basemap
bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)
# import shape file
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)

plt.show()


# fuelmask = np.zeros(len(WGSx))
# for nPoly in bm.fire_bounds_info:
# 	poly = polygons[nPoly]
# 	poly = zip(*poly)
# 	poly_path = path.Path(poly)
# 	mask = poly_path.contains_points(dem_lcn)
# 	if poly_type[nPoly]==1:
# 		fuelmask[mask==1] = 1
# 	elif poly_type[nPoly]==2:
# 		fuelmask[mask==1] = 0
# dem_landmask = np.reshape(fuelmask, np.shape(lon_grid))




# # apply to grid