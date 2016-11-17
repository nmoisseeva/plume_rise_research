# evaluation of fire spread and heat flux


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
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
wrfdata = '/Users/nadya2/data/plume/RxCADRE/wrfout_LG2'
wrfinput='/Users/nadya2/Applications/WRFV3/test/em_fire/wrfinput_d01'
bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'
instruments_shape = '/Users/nadya2/data/RxCADRE/instruments/HIP1'

target_ros = {'HIP1':0.225, 'HIP2':0.443,'HIP3':0.233} 	#rates of spread from Butler2016 for L2G
# HIP1_locs = '/Users/nadya2/data/RxCADRE/instruments/L2G_HIP1.csv'

ll_utm = np.array([521620,3376766]) 	#lower left corner of the domain in utm
basemap_path = '/Users/nadya2/code/plume/RxCADRE/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])
#=================end of input===============

print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
nc_inputdata = netcdf.netcdf_file(wrfinput, mode ='r')

#create a UTM grid
FUTMx = nc_inputdata.variables['FXLONG'][0,:,:] + ll_utm[0]
FUTMy = nc_inputdata.variables['FXLAT'][0,:,:] + ll_utm[1]
UTMx = nc_inputdata.variables['XLONG'][0,:,:] + ll_utm[0]
UTMy = nc_inputdata.variables['XLAT'][0,:,:] + ll_utm[1]

#convert coordinate systems to something basemaps can read
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())
WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))
FWGSx, FWGSy= pyproj.transform(epsg26916,wgs84,FUTMx.ravel(),FUTMy.ravel())
FWLONG, FWLAT = np.reshape(FWGSx, np.shape(FUTMx)), np.reshape(FWGSy, np.shape(FUTMy))

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

#read shapefiles
plt.figure()
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)		
instruments = bm.readshapefile(instruments_shape, name='hip1')
hip1_lcn = []
for nPt,item in enumerate(bm.hip1_info):
	if item['Instrument'] == 'FireBehaviorPackage':
		hip1_lcn.append(bm.hip1[nPt])
	if item['Sensor_ID'] == 'FB22': 		#get index of fb22 for individual comparison 
		fb22_idx = len(hip1_lcn)-1
hip1_lcn = np.array(hip1_lcn)
# bm.scatter(hip1_lcn[:,0],hip1_lcn[:,1])
# plt.show()
plt.close()

#construct KDtree from idealized grid for fire and atmospheric mesh
grid_coord_fire = zip(FWGSy,FWGSx)
FgridTree = KDTree(grid_coord_fire)
Fdist, Fgrid_id = FgridTree.query(hip1_lcn[:,::-1]) 	#reorder columnets to lat/lon
grid_coord = zip(WGSy,WGSx)
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(hip1_lcn[:,::-1]) 	#reorder columnets to lat/lon

#calculate average rate of spread
ros = np.copy(nc_data.variables['ROS'][:,:,:])
rosnan = ros	
rosnan[rosnan==0] = np.nan 			#mask all non-fire cells
rosnan = np.nanmean(rosnan,0)		#git time averaged values
l2g_ros = np.nanmean(np.nanmean(rosnan,0)) #get average value for the entire fire
print('Average ROS within fire area: %.2f m/s' %l2g_ros)

#calculate average peak heat flux
hfx = np.copy(nc_data.variables['GRNHFX'][:,:,:]) 	#extract fire heat flux
hfxnan = hfx 	
hfxnan[hfxnan<5] = np.nan 			#residence time defined as at least 5kW/m2 as per Butler2013
hfxnanmax = np.nanmax(hfxnan,0)/1000 	#get peak value in kW/m2
l2g_hfx_max = np.nanmean(np.nanmean(hfxnanmax,1)) #get average value for the entire fire
print('Average peak heat flux of the fire: %.2f kW m-2' %l2g_hfx_max)

#calculate average heat flux
hfxnanmean = np.nanmean(hfxnan,0)/1000 	#get peak value in kW/m2
l2g_hfx = np.nanmean(np.nanmean(hfxnanmean,1)) #get average value for the entire fire
print('Average heat flux of the fire: %.2f kW m-2' %l2g_hfx)


#calculate select values for HIP1
print('Calculating values for HIP1...')
ros_val, hfx_val, hfxmax_val = [],[],[]
for pt in Fgrid_id:
	flat_ros = rosnan.ravel() 
	pt_ros = flat_ros[pt]
	ros_val.append(pt_ros) 
	print('---> Probe val: ROS = %s m/s' %(pt_ros))
for nPt, pt in enumerate(grid_id):
	flat_hfx = hfxnanmean.ravel()
	flat_hfxmax = hfxnanmax.ravel()
	pt_hfx = flat_hfx[pt]
	pt_hfxmax = flat_hfxmax[pt]
	hfx_val.append(pt_hfx)
	hfxmax_val.append(pt_hfxmax)
	print('---> Probe val:  HFX = %s kW/m2, HFXmax = %s kW/m2' %(pt_hfx,pt_hfxmax))
	if nPt == fb22_idx:
		fb22_coords = np.unravel_index(pt, np.shape(hfx[0,:,:]))
		print('------> FB22 sensor: 3d grid IDs: y=%s x=%s' %fb22_coords )
mean_hip1_ros = np.mean(ros_val)
mean_hip1_hfx = np.mean(hfx_val)
mean_hip1_hfxmax = np.mean(hfxmax_val)
print('Averages for HIP1: ROS = %s m/s, HFX = %s kW/m2, HFXmax = %s kW/m2' %(mean_hip1_ros,mean_hip1_hfx,mean_hip1_hfxmax))


# ------------------------------PLOTTING-----------------------------
#plot mean ROS for the fire
plt.figure()
im = plt.contourf(rosnan) 	
plt.colorbar()
plt.title('ROS [m/s]')
plt.show()

#plot peak heat flux
plt.figure()
im = plt.contourf(hfxnanmax) 
plt.title('PEAK HFX DURING FIRE')
plt.colorbar()			
plt.show()

#plot average heat flux
plt.figure()
im = plt.contourf(hfxnanmean) 			
plt.colorbar()
plt.title('AVE HFX DURING FIRE')
plt.show()

#plot heat flux over FB22 sensor throughout the fire
plt.figure()
plt.plot(nc_data.variables['XTIME'][:],hfx[:,fb22_coords[0],fb22_coords[1]])
plt.title('HEAT FLUX AT FB22 SENSOR')
plt.xlabel('time [min]')
plt.ylabel('heat flux [W/m2]')
plt.show()






