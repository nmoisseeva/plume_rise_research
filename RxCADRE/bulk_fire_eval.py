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

sr_xy = 10								#fire grid downscale

ll_utm = np.array([519500,3377000])		#lower left corner of the domain in utm

basemap_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  

WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]
ghfx = np.copy(nc_data.variables['FGRNHFX'][:,:,:]) 	#extract fire heat flux
tsec = nc_data.variables['XTIME'][:] * 60 		#get time in seconds since run start
ros = np.copy(nc_data.variables['ROS'][:,:,:])


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
gdal.UseExceptions()

for nIm in range(num_im):
	#get closest time snapshot
	tend = mosaic_stamps[nIm,2]
	tend_sec = (int(tend[3:5]) - 27)*60 + int(tend[6:8])
	ti = np.argmin(abs(tsec - tend_sec))
	t_idx.append(ti)
 
	#get model heat flux
	hxsnap = ghfx[ti,:,:]
	hxsnap[hxsnap<5000] = np.nan

	# plt.contourf(hxsnap)
	# plt.show()
	# plt.hist(hxsnap[np.isfinite(hxsnap)], bins=20)
	# plt.show()

	#load tiff mosaic files
	trange = mosaic_stamps[nIm,0]
	print('Current mosaic: %s' %trange)


	outtif = lwir_data + 'NMcrop_' + trange + '.tif'
	if not os.path.isfile(outtif):
		intif = lwir_data + 'L2G_Wm2_mosaic_' + trange + '.tif'
		warp_cmd = 'gdalwarp -cutline ' +  bounds_shape + ' -crop_to_cutline -dstalpha ' +  intif + ' ' + outtif
		os.system(warp_cmd)
	lwir = gdal.Open(outtif)
	band = lwir.GetRasterBand(1)
	tiffdata = np.array(band.ReadAsArray(),dtype = float)
	tiffdata[tiffdata==65535] = np.nan
	tiffdata[tiffdata<5000] = np.nan
	# plt.imshow(tiffdata)
	# plt.colorbar()
	# plt.show()

	#get resolution of LIWR and model 
	geotransform = lwir.GetGeoTransform()					
	pixTif = abs(geotransform[1]) * abs(geotransform[5])
	pixMod = nc_data.DX * nc_data.DY / (sr_xy**2)

	#get mean fluxes
	aveHfxTif, aveHfxMod = np.nanmean(tiffdata.ravel()), np.nanmean(hxsnap.ravel())

	#get ros values and averages
	rossnap = ros[ti,:,:]
	rossnap[np.isnan(hxsnap)] = np.nan
	aveRosMod = np.nanmean(rossnap.ravel())
	print('.....Average model ROS: %.2f (m/s2)' %aveRosMod)

	#get total fire output
	flat_lwir_vals = tiffdata[~np.isnan(tiffdata)] * pixTif
	flat_mod_vals = hxsnap[~np.isnan(hxsnap)] * pixMod
	areaTif = pixTif * len(flat_lwir_vals)
	areaMod = pixMod * len(flat_mod_vals)
	totHeatTif, totHeatMod = np.sum(flat_lwir_vals)/1000000., np.sum(flat_mod_vals)/1000000.

	print('.....Avearage heat flux real vs model fire: %.2f vs. %.2f (W/m2)' %(aveHfxTif, aveHfxMod))
	print('.....Total heat output of real vs model fire: %.2f vs. %.2f (MW)' %( totHeatTif,totHeatMod))
	print('.....Total fire area of real vs model fire: %.2f vs. %.2f (ha)' %( areaTif*1e-4,areaMod*1e-4))
