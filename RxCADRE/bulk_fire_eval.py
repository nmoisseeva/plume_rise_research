import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
from scipy.spatial import KDTree
from matplotlib import path 
from mpl_toolkits import basemap
import mpl_toolkits.basemap.pyproj as pyproj
import os.path
import pickle
import gdal


#====================INPUT===================
wrfdata = '/Users/nmoisseeva/data/plume/RxCADRE/regrid/wrfout_L2G_cat3_new_regrid'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
bounds_shape = '/Users/nmoisseeva/data/qgis/active_burn_perim.shp'
lwir_data = '/Users/nmoisseeva/data/RxCADRE/LWIR/RDS-2016-0008/Data/2012_L2G/'
mosaic_data = '/Users/nmoisseeva/data/RxCADRE/LWIR/MosaicTimes_NM.csv'


ll_utm = np.array([519500,3377000])		#lower left corner of the domain in utm

basemap_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  

WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]
ghfx = np.copy(nc_data.variables['GRNHFX'][:,:,:]) 	#extract fire heat flux
tsec = nc_data.variables['XTIME'][:] * 60 		#get time in seconds since run start


#open/generate basemap
if os.path.isfile(basemap_path):
	bm = pickle.load(open(basemap_path,'rb'))   # load here the above pickle 
	print('Domain basemap found at: %s' %basemap_path)
else:
	print('WARNING: no existing basemaps found: configuring a new basemap')
	bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)
	pickle.dump(bm,open(basemap_path,'wb'),-1)  	# pickle the new map for later
	print('.....New basemap instance saved as: %s' %basemap_path)


#load mosaic start and end times
mosaic_stamps = np.genfromtxt(mosaic_data, skip_header=1, delimiter=',',dtype=str)
num_im = np.shape(mosaic_stamps)[0]

t_idx = []
for nIm in range(num_im):
	tend = mosaic_stamps[nIm,2]
	tend_sec = (int(tend[3:5]) - 27)*60 + int(tend[6:8])
	ti = np.argmin(abs(tsec - tend_sec))
	t_idx.append(ti)
	hxsnap = ghfx[ti,:,:]
	hxsnap[hxsnap==0] = np.nan
	plt.contourf(hxsnap)
	plt.show()
	print np.nanmean(hxsnap.ravel())

gdal.UseExceptions()
for nIm in range(num_im):
	#load tiff mosaic files starting from the last snapshot
	trange = mosaic_stamps[nIm,0]
	intif = lwir_data + 'L2G_Wm2_mosaic_' + trange + '.tif'
	outtif = lwir_data + 'NMcrop_' + trange + '.tif'
	if not os.path.isfile(outtif):
		warp_cmd = 'gdalwarp -cutline ' +  bounds_shape + ' -crop_to_cutline -dstalpha ' +  intif + ' ' + outtif
		os.system(warp_cmd)
	lwir = gdal.Open(outtif)
	band = lwir.GetRasterBand(1)
	tiffdata = np.array(band.ReadAsArray(),dtype = float)
	tiffdata[tiffdata==65535] = np.nan
	tiffdata[tiffdata<600] = np.nan
	plt.imshow(tiffdata)
	plt.colorbar()
	plt.show()
	print np.nanmean(tiffdata.ravel())




# gdal.UseExceptions()
# for nIm in range(num_im):
# 	#load tiff mosaic files starting from the last snapshot
# 	trange = mosaic_stamps[nIm,0]
# 	geotiff = lwir_data + 'L2G_Wm2_mosaic_' + trange + '.tif'
# 	lwir = gdal.Open(geotiff)
# 	band = lwir.GetRasterBand(1)
# 	tiffdata = band.ReadAsArray()

# 	#create lat/lon grids
# 	gt = lwir.GetGeoTransform()							#get raster data			
# 	nPy, nPx = np.shape(tiffdata)						#get size of data
# 	lat_grid = np.zeros((nPx,nPy))						#create storage arrays
# 	lon_grid = np.zeros((nPx,nPy))
# 	for nY in range(nPy):
# 		for nX in range(nPx):
# 			lat_grid[nX,nY] = gt[3] + nX*gt[4] + nY*gt[5]
# 			lon_grid[nX,nY] = gt[0] + nX*gt[1] + nY*gt[2]

# 	min_lat_idx = np.argmin(np.abs(lat_grid[0,:]-lat[0])) 
# 	max_lat_idx = np.argmin(np.abs(lat_grid[0,:]-lat[1]))
# 	min_lon_idx = np.argmin(np.abs(lon_grid[:,0]-lon[0]))
# 	max_lon_idx = np.argmin(np.abs(lon_grid[:,0]-lon[1]))

# 	bounds = (lon_grid[min_lon_idx,0],lon_grid[max_lon_idx,0], lat_grid[0,min_lat_idx],lat_grid[0,max_lat_idx])

# 	if any((min_lat_idx,max_lat_idx)) == any((0,nPy)):
# 		print("WARNING: requestested latitude range extends to the edge (or beyond) the available geoTIFF domain")
# 	if any((min_lat_idx,max_lat_idx)) == any((0,nPx)):
# 		print("WARNING: requestested longitude range extends to the edge (or beyond) the available geoTIFF domain")

# 	minY, maxY = min(max_lat_idx,min_lat_idx), max(max_lat_idx,min_lat_idx)
# 	minX, maxX = min(max_lon_idx,min_lon_idx), max(max_lon_idx,min_lon_idx)
# 	lat_grid = lat_grid[minX:maxX, minY:maxY]
# 	lon_grid = lon_grid[minX:maxX, minY:maxY]
# 	tiffdata = tiffdata[minY:maxY, minX:maxX]

# 	tiffdata  = tiffdata.T


