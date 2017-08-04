# evaluation of fire spread and heat flux


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
from scipy.ndimage.interpolation import zoom
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
# wrfdata = '/Users/nadya2/data/plume/RxCADRE/regrid/wrfout_LG2_nospinup_regrid'
wrfdata = '/Users/nmoisseeva/data/plume/RxCADRE/regrid/wrfout_L2G_cat1_regrid'

bounds_shape = '/Users/nmoisseeva/data/qgis/LG2012_WGS'
instruments_shape = '/Users/nmoisseeva/data/RxCADRE/instruments/HIP1'

# burn_lmt = [(-86.73051,30.54315),(-86.73398,30.54584),(-86.74677,30.53858),(-86.74422,30.53546)]

ign_lcn = (524600,3378360)

# ll_utm = np.array([518800,3377000])		#lower left corner of the domain in utm
ll_utm = np.array([519500,3377000])		#lower left corner of the domain in utm

basemap_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])
#=================end of input===============

print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  

tstep = round((nc_data.variables['XTIME'][1]- nc_data.variables['XTIME'][0]) * 60.)  #timestep in sec
WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]

#convert coordinate systems to something basemaps can read
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")


FWGSx, FWGSy = zoom(WLONG,10),zoom(WLAT,10)
FUTMx, FUTMy = pyproj.transform(wgs84,epsg26916,FWGSx.ravel(),FWGSy.ravel())

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
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)		
instruments = bm.readshapefile(instruments_shape, name='hip1')
hip1_lcn, hip1_tag = [], []
for nPt,item in enumerate(bm.hip1_info):
	if item['Instrument'] == 'FireBehaviorPackage':
		hip1_lcn.append(bm.hip1[nPt])
		hip1_tag.append(item['Sensor_ID'])

hip1_lcn, hip1_tag = np.array(hip1_lcn), np.array(hip1_tag)
# bm.contourf(FWGSx,FWGSy,fhfx[100,:,:])
plt.show()

#convert hip1 locations to UTM
hipUTMx, hipUTMy = pyproj.transform(wgs84,epsg26916, hip1_lcn[:,0],hip1_lcn[:,1])

print('Create a KDtree of fire')
grid_coord_fire = zip(FUTMy,FUTMx)
FgridTree = KDTree(grid_coord_fire)
Fdist, Fgrid_id = FgridTree.query(zip(hipUTMy, hipUTMx)) 	
Idist, Igrid_id = FgridTree.query([ign_lcn[1], ign_lcn[0]]) 	#location of ignition line normal

#get necessary variables
ghfx = np.copy(nc_data.variables['GRNHFX'][:,:,:]) 	#extract fire heat flux
fhfx = np.copy(nc_data.variables['FGRNHFX'][:,:-10,:-10]) 	#extract fire heat flux
fuelfrac = np.copy(nc_data.variables['AVG_FUEL_FRAC'][:,:,:])
xtime = nc_data.variables['XTIME'][:] * 60 			#get time in seconds (since noon)

#create ignition mask (fire mesh)
print('..... creating an igntion mask (may take several minutes)')
ign_mask = np.empty_like(fhfx) * np.nan
for nt in range(len(xtime)):
	print nt
	current_ign = fhfx[nt,:,:]
	temp_mask = np.empty_like(current_ign) * np.nan
	temp_mask[current_ign>5000] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	# temp_mask[current_ign>0] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013

	ign_mask[nt,:,:] = temp_mask

print('..... creating an igntion mask on atm grid (may take several minutes)')
ign_mask_atm = np.empty_like(ghfx) * np.nan
for nt in range(len(xtime)):
	print nt
	frac = fuelfrac[nt,:,:]
	current_ign = ghfx[nt,:,:]
	temp_mask = np.empty_like(current_ign) * np.nan
	temp_mask[current_ign>5000] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	# temp_mask[current_ign>0] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	# temp_mask[(frac<1) & (frac>0)] = 1 
	ign_mask_atm[nt,:,:] = temp_mask

#calculate average peak heat flux on fire grid
ignFHFX = np.copy(fhfx) 			#convert to kW/m2
ignFHFX[np.isnan(ign_mask)] = np.nan 	#get ignited cells only

aveFHFX = np.nanmean(ignFHFX,0)/1000. #get ave value in kW/m2
maxFHFX = np.nanmax(ignFHFX,0)/1000. #get peak value in kW/m2

#calculate average peak heat flux on atmospheric grid
ignGHFX = np.copy(ghfx) 			#convert to kW/m2
ignGHFX[np.isnan(ign_mask_atm)] = np.nan 	#get ignited cells only

aveGHFX = np.nanmean(ignGHFX,0)/1000. #get ave value in kW/m2
maxGHFX = np.nanmax(ignGHFX,0)/1000. #get peak value in kW/m2

print('Average average heat flux of the fire: %.2f kW m-2' %np.nanmean(aveGHFX))
print('Average peak flux of the fire: %.2f kW m-2' %np.nanmean(maxGHFX))

#calculate values for HIP1
print('Calculating values for HIP1...')
hip1_max, hip1_ave, t_ign = [], [], []
hip1_hfx = np.empty((len(hip1_tag),len(xtime)))*np.nan

for nPt, pt in enumerate(Fgrid_id):
	#get max values for sensors
	flat_max = maxFHFX.ravel() 
	pt_hfxmax = flat_max[pt]
	hip1_max.append(pt_hfxmax) 

	#get average values for sensors
	pt_idx = np.unravel_index(pt,np.shape(FWGSx))
	print pt_idx
	hip1_hfx[nPt,:] = ignFHFX[:,pt_idx[0],pt_idx[1]]/1000.
	pt_ave = np.nanmean(hip1_hfx[nPt,:])
	hip1_ave.append(pt_ave)
	print('---> Sensor %s: hfx = %s kW/m2, maxhfx = %s kW/m2' %(hip1_tag[nPt],pt_ave,pt_hfxmax))

	#calculate rate of spread
	pt_ign = np.argmax(np.isfinite(hip1_hfx[nPt,:]))  	#get index of first non-nan
	t_ign.append(pt_ign)

ign_idx = np.unravel_index(Igrid_id,np.shape(FWGSx)) 				#get location of normal
t0 = np.argmax(np.isfinite(ignFHFX[:,ign_idx[0],ign_idx[1]])) 	#get time of igntion of normal
t = (t_ign - t0) * tstep
print('.....fireline igntion t0 frame: %s' %t0)
print('.....ignition time in sec: %s' %t_ign)

subTree = KDTree(zip(hipUTMx,hipUTMy))
rosdist, ros_id = subTree.query((ign_lcn[0], ign_lcn[1]),k=len(hip1_tag)) 	#reorder columnets to lat/lon

ros = rosdist / t
aveROS = np.mean(ros)
print('.....Avarage rate of spread based on HIP1: %s' %aveROS)

#--------------------------DATA from Butler 2013--------------------
butler_data = {'aveHFX':{}, 'maxHFX':{}, 'ROS':{}}

butler_data['ROS']['arrival'] = [0.32,0.24,0.25,0.22,0.22,0.10]
butler_data['ROS']['FRP'] = [0.24,0.31,0.44,0.61]

butler_data['maxHFX']['HIP1'] = [36.7,32.3,12.5,30.8,13.6,8.8,0.7]
butler_data['maxHFX']['L2G'] = [36.7,32.3,12.5,30.8,13.6,8.8,0.7,12.9,26.9,23.5,25.5,22.1,24,13.2,20,30.2,13.1,16.2,29.8,17.7,8.6,]

butler_data['aveHFX']['HIP1'] = [19.3,16.9,5.9,14.1,0.7]
butler_data['aveHFX']['L2G'] = [19.3,16.9,5.9,14.1,0.7,5.2,15,13.9,15,14.2,15.6,8.7,8.6,12.8,6.2,8.6,13.3,8.3,6.3]

# ------------------------------PLOTTING-----------------------------
plt.figure()
plt.title('HEAT FLUX AT INDIVIDUAL HIP1 SENSORS')
for i in range(len(hip1_tag)):
	plt.plot(hip1_hfx[i,:], label=hip1_tag[i])
plt.xticks(np.arange(0,len(xtime),20),xtime[::20])
plt.xlabel('time [sec]')
plt.ylabel('heat flux $[kW m^{-2}]$')
plt.legend()
plt.show()

plt.figure()
plt.title('AVERAGE HEAT FLUX')
box = plt.boxplot([butler_data['aveHFX']['L2G'],\
					butler_data['aveHFX']['HIP1'],\
					aveGHFX[np.isfinite(aveGHFX)],\
					hip1_ave], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES','LES HIP1'])
plt.ylabel('heat flux $[kW m^{-2}]$')
plt.show()

plt.figure()
plt.title('PEAK HEAT FLUX')
box = plt.boxplot([butler_data['maxHFX']['L2G'],\
					butler_data['maxHFX']['HIP1'],\
					aveGHFX[np.isfinite(maxGHFX)],\
					hip1_max], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES','LES HIP1'])
plt.ylabel('heat flux $[kW m^{-2}]$')
plt.show()


plt.figure()
plt.title('ROS')
box = plt.boxplot([butler_data['ROS']['arrival'],\
					butler_data['ROS']['FRP'],\
					ros], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES'])
plt.ylabel('ROS $[m s^{-1}]$')
plt.show()

# #plot mean ROS for the fire
# plt.figure()
# im = plt.contourf(rosnan[0:1000,1000:]) 	
# plt.colorbar()
# plt.title('ROS [m/s]')
# plt.show()

# #plot peak heat flux
# plt.figure()`
# im = plt.contourf(hfxnanmax[0:100,100:]) 
# plt.title('PEAK HFX DURING FIRE')
# plt.colorbar()			
# plt.show()

# #plot average heat flux
# plt.figure()
# im = plt.contourf(hfxnanmean[0:100,100:]) 			
# plt.colorbar()
# plt.title('AVE HFX DURING FIRE')
# plt.show()




