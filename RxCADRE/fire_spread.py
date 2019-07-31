# evaluation of fire spread and heat flux


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy.spatial import KDTree
from scipy.ndimage.interpolation import zoom
import matplotlib.animation as animation
from matplotlib import path
from mpl_toolkits import basemap
# import mpl_toolkits.basemap.pyproj as pyproj
import shapefile
import pyproj
import os.path
import pickle
import mpl_toolkits.mplot3d.axes3d as p3
import mpl_toolkits.mplot3d as a3
from matplotlib import animation
import imp

#====================INPUT===================
#all common variables are stored separately
import rxcadreMOIST as rx
imp.reload(rx)		        #force load each time

#==============calculated input===============
#
nsamples = 20
#ignition start samples above 2nd and 3rd lines
ign_lcn_l2g_x = np.append(np.linspace(525527.,524584.,nsamples),np.linspace(525432.,524517.,nsamples))
ign_lcn_l2g_y = np.append(np.linspace(3378993.,3378385.,nsamples),np.linspace(3379091.,3378503.,nsamples))
#end samples above 2nd and 3rd lines
test_lcn_l2g_x = np.append(np.linspace(525474.,524540.,nsamples),np.linspace(525390.,524475.,nsamples))
test_lcn_l2g_y = np.append(np.linspace(3379050.,3378430.,nsamples),np.linspace(3379135.,3378547.,nsamples))


#define normals to fire spread for HIP1 - in order ['FB19', 'FB2', 'FB3', 'FB14', 'FB17', 'FB20', 'FB22']
ign_lcn_hip_x = [524604,524621,524607,524606,524604,524594,524577]
ign_lcn_hip_y = [3378359,3378377,3378366,3378365,3378360,3378352,3378339]

#original ignition line locations (from namelist.input)
wrf_lines_x_start = [8507.,8630.,8732.,8820.] + rx.ll_utm[0]
wrf_lines_x_end = [7305.,7384.,7462.,7548.] + rx.ll_utm[0]
wrf_lines_y_start = [2248.,2179.,2077.,2004.] + rx.ll_utm[1]
wrf_lines_y_end = [1463.,1366.,1248.,1172.] + rx.ll_utm[1]

#=================end of input===============
print('FIRE BEHAVIOR ANALYSIS')

print('.....extracting NetCDF data from %s ' %rx.wrfdata)
ncdata = netcdf.netcdf_file(rx.wrfdata, mode ='r')

#load georeferencing data for the same run
print('.....importing wrf coordinates from  %s ' %rx.geo_path)
wrfgeo = np.load(rx.geo_path, allow_pickle=True).item()

tstep = round((ncdata.variables['XTIME'][1]- ncdata.variables['XTIME'][0]) * 60.)  #timestep in sec
# WLONG, WLAT = ncdata.variables['XLONG'][0,:,:], ncdata.variables['XLAT'][0,:,:]
# #
# #convert coordinate systems to something basemaps can read
# wgs84=pyproj.Proj("+init=EPSG:4326")
# epsg26916=pyproj.Proj("+init=EPSG:26916")
#
# FWGSx, FWGSy = zoom(WLONG,10),zoom(WLAT,10)
# FUTMx, FUTMy = pyproj.transform(wgs84,epsg26916,FWGSx.ravel(),FWGSy.ravel())
# UTMx, UTMy = pyproj.transform(wgs84,epsg26916,WLONG.ravel(),WLAT.ravel())

#open/generate basemap
if os.path.isfile(rx.basemap_path):
	bm = pickle.load(open(rx.basemap_path,'rb'))   # load here the above pickle
	print('Domain basemap found at: %s' %rx.basemap_path)
else:
	print('WARNING: no existing basemaps found: configuring a new basemap')
	bm = basemap.Basemap(llcrnrlon=WLONG[0,0], llcrnrlat=WLAT[0,0],\
					 urcrnrlon=WLONG[-1,-1], urcrnrlat=WLAT[-1,-1], resolution='f', epsg=4326)
	pickle.dump(bm,open(rx.basemap_path,'wb'),-1)  	# pickle the new map for later
	print('.....New basemap instance saved as: %s' %rx.basemap_path)

#read shapefiles
polygons = bm.readshapefile(rx.bounds_shape,name='fire_bounds',drawbounds=True)
instruments = bm.readshapefile(rx.hip1_shape, 'hip1')
hip1_tag, hip1UTMx, hip1UTMy = [], [], []
for nPt,item in enumerate(bm.hip1_info):
    if item['Instrument'] == 'FireBehaviorPackage':
        hip1_tag.append(item['Sensor_ID'])
        hip1UTMy.append(item['Northing_Y'])
        hip1UTMx.append(item['Easting_X'])

print('..... creating a KDtree of fire mesh')
grid_coord_fire = list(zip(wrfgeo['FXLAT'].ravel(),wrfgeo['FXLONG'].ravel()))
FgridTree = KDTree(grid_coord_fire)
Fdist, Fgrid_id = FgridTree.query(list(zip(hip1UTMy, hip1UTMx)))
Idist, Igrid_id = FgridTree.query(list(zip(ign_lcn_hip_y,ign_lcn_hip_x))) 	#location of ignition line normal
TESTdist, TESTgrid_id = FgridTree.query(list(zip(test_lcn_l2g_y,test_lcn_l2g_x)))
NORMdist, NORMgrid_id = FgridTree.query(list(zip(ign_lcn_l2g_y,ign_lcn_l2g_x)))

print('..... creating a KDtree of atm mesh')
grid_coord_atm = list(zip(wrfgeo['XLAT'].ravel(),wrfgeo['XLONG'].ravel()))
gridTree = KDTree(grid_coord_atm)
Adist, Agrid_id = gridTree.query(list(zip(hip1UTMy, hip1UTMx)))

#get necessary variables
ghfx = ncdata.variables['GRNHFX'][:,:,:] 	          #extract fire heat flux
fhfx = ncdata.variables['FGRNHFX'][:,:,:] 	#extract fire heat flux
xtime = ncdata.variables['XTIME'][:] * 60 			#get time in seconds

#create ignition mask (fire mesh)
print('..... creating an igntion mask (may take several minutes)')
ign_mask = np.empty_like(fhfx) * np.nan
for nt in range(len(xtime)):
	print(nt)
	current_ign = fhfx[nt,:,:]
	temp_mask = np.empty_like(current_ign) * np.nan
	temp_mask[current_ign>5000] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	ign_mask[nt,:,:] = temp_mask

print('..... creating an igntion mask on atm grid (may take several minutes)')
ign_mask_atm = np.empty_like(ghfx) * np.nan
for nt in range(len(xtime)):
	print(nt)
	current_ign = ghfx[nt,:,:]
	temp_mask = np.empty_like(current_ign) * np.nan
	temp_mask[current_ign>5000] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	ign_mask_atm[nt,:,:] = temp_mask

print("Calculating peak and averages heat flux values for all ignited cells:")
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

print('..... Average average heat flux of the fire (all cells): %.2f kW m-2' %np.nanmean(aveGHFX))
print('..... Average peak flux of the fire (all cells): %.2f kW m-2' %np.nanmean(maxGHFX))

print("----------------HIP1 analysis------------------")
print("Calculating HIP1 sensor heat fluxes on fire mesh:")
#calculate values for HIP1
hip1_max, hip1_ave, t_ign_hip1 = [], [], []
hip1_hfx = np.empty((len(hip1_tag),len(xtime)))*np.nan
#----on fire mesh (currentlly only used for ROS)
for nPt, pt in enumerate(Fgrid_id):
	#get max values for sensors
	flat_max = maxFHFX.ravel()
	pt_hfxmax = flat_max[pt]
	hip1_max.append(pt_hfxmax)

	#get average values for sensors
	pt_idx = np.unravel_index(pt,np.shape(wrfgeo['FXLONG']))
	hip1_hfx[nPt,:] = ignFHFX[:,pt_idx[0],pt_idx[1]]/1000.
	pt_ave = np.nanmean(hip1_hfx[nPt,:])
	hip1_ave.append(pt_ave)
	print('..... Sensor %s: avehfx = %s kW/m2, maxhfx = %s kW/m2' %(hip1_tag[nPt],pt_ave,pt_hfxmax))

	#for calculating rate of spread
	pt_ign = np.argmax(np.isfinite(hip1_hfx[nPt,:]))  	#get index of first non-nan
	t_ign_hip1.append(pt_ign)

print("Calculating HIP1 sensor heat fluxes on atmospheric mesh:")
#----on atm grid (used for heat fluxes)
hip1_max_atm, hip1_ave_atm  = [],[]
for nPt, pt in enumerate(Agrid_id):
	#get max values for sensors
	flat_max = maxGHFX.ravel()
	pt_hfxmax = flat_max[pt]
	hip1_max_atm.append(pt_hfxmax)

	#get average values for sensors
	pt_idx = np.unravel_index(pt,np.shape(wrfgeo['WLAT']))
	hip1_hfx[nPt,:] = ignGHFX[:,pt_idx[0],pt_idx[1]]/1000.
	pt_ave = np.nanmean(hip1_hfx[nPt,:])
	hip1_ave_atm.append(pt_ave)
	print('..... Sensor %s: hfx = %s kW/m2, maxhfx = %s kW/m2' %(hip1_tag[nPt],pt_ave,pt_hfxmax))

print('..... Average average heat flux of the fire (HIP1): %.2f kW m-2' %np.nanmean(hip1_ave_atm))
print('..... Average peak flux of the fire (HIP1): %.2f kW m-2' %np.nanmean(hip1_max_atm))

print("Calculating ROS for HIP1:")
# #calculate ROS for hip1
dt_hip = []
for nPt, pt in enumerate(Igrid_id):
	ign_idx = np.unravel_index(pt,np.shape(wrfgeo['FXLONG'])) 		#get index of normal location
	t0 = np.argmax(np.isfinite(ignFHFX[:,ign_idx[0],ign_idx[1]]))	#get time of igntion of normal
	dt_hip.append((t_ign_hip1[nPt]-t0) * tstep)
dist_hip = np.sqrt((ign_lcn_hip_x - np.array(hip1UTMx))**2 + (ign_lcn_hip_y - np.array(hip1UTMy))**2)
ros_hip = dist_hip/dt_hip
aveROS_hip = np.mean(ros_hip)
print('.....Avarage rate of spread based on HIP1: %s' %aveROS_hip)


print("----------------Mid-section sample analysis------------------")
print("Calculating ROS based on cross sections:")
test_lcn_hfx = np.empty((len(ign_lcn_l2g_y),len(xtime))) * np.nan
t_ign_l2g = [] 				#time of sample ignition storage
for nPt, pt in enumerate(TESTgrid_id):
	pt_idx = np.unravel_index(pt,np.shape(wrfgeo['FXLONG'])) 				#get index of sample location
	test_lcn_hfx[nPt,:] = ignFHFX[:,pt_idx[0],pt_idx[1]]/1000	#get heat flux at sample
	if any(np.isfinite(test_lcn_hfx[nPt,:])):
		pt_ign = np.where(np.isfinite(test_lcn_hfx[nPt,:]))[0][0]  	#get index of first non-nan
	else:
		pt_ign = np.nan
	t_ign_l2g.append(pt_ign)

#calculate ROS along midfire
dt_l2g = []
dist_l2g = np.sqrt((ign_lcn_l2g_x - test_lcn_l2g_x)**2 + (ign_lcn_l2g_y - test_lcn_l2g_y)**2)
for nPt, pt in enumerate(NORMgrid_id):
	norm_idx = np.unravel_index(pt,np.shape(wrfgeo['FXLONG']))
	t0 = np.argmax(np.isfinite(ignFHFX[:,norm_idx[0],norm_idx[1]])) #get time of igntion of normal
	dt_pt = (t_ign_l2g[nPt] - t0) * tstep
	dt_l2g.append(dt_pt)
ros_l2g = dist_l2g/dt_l2g
aveROS_l2g = np.nanmean(ros_l2g)
print('.....Avarage rate of spread based on cross-section: %s' %aveROS_l2g)

#--------------------------DATA from Butler 2013--------------------
butler_data = {'aveHFX':{}, 'maxHFX':{}, 'ROS':{}}

butler_data['ROS']['L2G'] = [0.32,0.24,0.25,0.22,0.22,0.10,0.33,0.46,0.47,0.47,0.44,0.45,0.48,0.25,0.39,0.16,0.30,0.27,0.16,0.10]
butler_data['ROS']['HIP1'] = [0.32,0.24,0.25,0.22,0.22,0.10]

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
					hip1_ave_atm], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES','LES HIP1'])
plt.ylabel('heat flux $[kW m^{-2}]$')
plt.savefig(rx.fig_dir + 'AveHx.pdf')
plt.show()
plt.close()

plt.figure()
plt.title('PEAK HEAT FLUX')
box = plt.boxplot([butler_data['maxHFX']['L2G'],\
					butler_data['maxHFX']['HIP1'],\
					maxGHFX[np.isfinite(maxGHFX)],\
					hip1_max_atm], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES','LES HIP1'])
plt.ylabel('heat flux $[kW m^{-2}]$')
plt.savefig(rx.fig_dir + 'MaxHx.pdf')
plt.show()
plt.close()


plt.figure()
plt.title('ROS')
box = plt.boxplot([butler_data['ROS']['L2G'],\
					butler_data['ROS']['HIP1'],\
					ros_l2g,\
					ros_hip], notch=True, patch_artist=True,showmeans=True)
colors = ['lightblue','lightblue','pink','pink']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.xticks([1,2,3,4],['L2G','HIP1','LES','LES HIP1'])
plt.ylabel('ROS $[m s^{-1}]$')
plt.ylim([0,0.5])
plt.savefig(rx.fig_dir + 'ROS.pdf')
plt.show()
plt.close()

plt.figure()
# x = np.reshape(FUTMx, np.shape(FWGSx))
# y = np.reshape(FUTMy, np.shape(FWGSy))
plt.contourf(wrfgeo['FXLONG'][250:600,-900:-300],wrfgeo['FXLAT'][250:600,-900:-300],fhfx[100,250:600,-900:-300],cmap=plt.cm.gist_heat_r)
plt.scatter(ign_lcn_l2g_x,ign_lcn_l2g_y, s=2, label='initial sample')
plt.scatter(test_lcn_l2g_x,test_lcn_l2g_y,s=2, label='final sample')
for nLine in range(4):
	plt.plot((wrf_lines_x_start[nLine],wrf_lines_x_end[nLine]),\
			(wrf_lines_y_start[nLine],wrf_lines_y_end[nLine]), \
			'r--',linewidth =1, label='WRF-SFIRE ignition lines' )
plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
from collections import OrderedDict
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.show()
