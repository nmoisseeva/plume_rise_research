# referencing simulation with flight data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
from scipy.ndimage.interpolation import rotate
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
wrfdata = '/Users/nadya2/data/plume/RxCADRE/regrid/wrfout_LG2'
fig_dir = '/Users/nadya2/code/plume/figs/RxCADRE/'
bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'
disp_data = '/Users/nadya2/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
emis_data = '/Users/nadya2/data/RxCADRE/dispersion/Data/Emissions_L2G_20121110.csv'

ll_utm = np.array([518800,3377000]) 	#lower left corner of the domain in utm
basemap_path = '/Users/nadya2/code/plume/RxCADRE/npy/%s_%s_bm.npy' %(ll_utm[0],ll_utm[1])

lvl = np.arange(0,2000,50) 				#
emis_excl = 0 							#number of samples to excluded (from the END!)
sfc_hgt = 62 							#surface height MSL (m)
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
# nc_interp = netcdf.netcdf_file(wrfinterp, mode ='r')  

#get geopotential array and convrt to height
z = (nc_data.variables['PHB'][:,:,:,:] + nc_data.variables['PH'][:,:,:,:]) / 9.81

# #create a UTM grid
# UTMx = nc_data.variables['XLONG'][0,:,:] + ll_utm[0]
# UTMy = nc_data.variables['XLAT'][0,:,:] + ll_utm[1]

# #convert coordinate systems to something basemaps can read
# wgs84=pyproj.Proj("+init=EPSG:4326")
# epsg26916=pyproj.Proj("+init=EPSG:26916")
# WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())
# WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))

WLONG, WLAT = nc_data.variables['XLONG'][0,:,:], nc_data.variables['XLAT'][0,:,:]


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

# Sanity check: import shape file
# polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)
# fireim = nc_data.variables['GRNHFX'][25,:,:]
# bm.imshow(fireim)
# plt.show()

#extract model time info
runstart = nc_data.START_DATE[-8:]
tsec = nc_data.variables['XTIME'][:] * 60 		#get time in seconds since run start
model_ssm = int(runstart[0:2])*3600 + int(runstart[3:5])*60

#==========================VERTICAL INTERPOLATION============================
qvcopy = np.copy(nc_data.variables['QVAPOR'][:,:,:,:])



#================================DISPERSION==================================
#extract and format dispersion data
disp_dict = {}
disp_array = np.genfromtxt(disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')

start_idx = np.argmin(abs(disp_array[:,0] - model_ssm))  #find simulation start time index
disp_dict['time']= disp_array[start_idx:,0] - model_ssm +1
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]-sfc_hgt))
disp_dict['meta']= 'time: seconds since model start run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'


#get indecies of samples corresponding to model output times
tidx = [np.argmin(abs(disp_dict['time']-t)) for t in tsec] #times since start
dt = disp_dict['time'][1] - disp_dict['time'][0] 	

# #construct KDtree from idealized grid
# grid_coord = zip(WGSy,WGSx,lvl)
# gridTree = KDTree(grid_coord)
# dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx])

#calculate H20 MR departure and extract point locations (assumes dry initial state)
# qvapor = np.copy(nc_interp.variables['QVAPOR'][:,:,:,:])
qvapor = np.copy(nc_data.variables['QVAPOR'][:,:,:,:])
qvapornan = np.copy(qvapor)
qvapornan[qvapornan<1e-30] = np.nan
mod_val = np.empty((len(tsec))) * np.nan
# obs_val = np.empty((len(tsec))) * np.nan
print('Finding nearest points....')
for nt in range(len(tsec)):
	print('tstep: %s') %nt
	#construct KDtree from idealized grid
	lvl = z[nt,:,:,:].ravel()
	grid_coord = zip(WGSy,WGSx,lvl)
	gridTree = KDTree(grid_coord)
	dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx])

	flatq = qvapornan[nt,::-1,:,:].ravel()
	mod_val[nt] = flatq[grid_id[nt]]
	mv = disp_dict['H2O'][tidx[nt]]		#in percent by volume
	# obs_val[nt] = 0.622*mv/(100-mv)

# plt.scatter(disp_dict['time'],disp_dict['CO2'])
plt.scatter(disp_array[:,0],disp_array[:,3])
# plt.plot(tsec,mod_val)
plt.show()

#================================EMISSIONS==================================
#extract and format emissions data 
emis_dict = {}
emis_array = np.genfromtxt(emis_data, skip_header=1, usecols = [5,6,13,14], delimiter=',',skip_footer=emis_excl)

emis_dict['bkgd']= zip(emis_array[:,0] - model_ssm +1, emis_array[:,1] - model_ssm +1)
emis_dict['smoke'] = zip((emis_array[:,2] - model_ssm +1).astype(int), (emis_array[:,3] - model_ssm +1).astype(int))
emis_dict['meta']= 'bkgd: background slice start and end in sec from simulation start | \
					smoke: plume start and end in sec from simulation start | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation MSL'



#================================PLOTTING==================================
#plot of CO2 slices
plt.figure()
plt.title('CO2 anomalies')
plt.scatter(disp_dict['time'],disp_dict['CO2'])
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 390,440, facecolor='gray', alpha=0.1, edgecolor='w')
plt.ylim([390,440])
plt.xlim([0,3000])
plt.show()


#plot of H2O slices
plt.figure()
plt.title('H2O anomalies')
plt.scatter(disp_dict['time'],disp_dict['H2O'])
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 0,1.5, facecolor='gray', alpha=0.1, edgecolor='w')
plt.ylim([0,1.5])
plt.xlim([0,3000])
plt.show()

#plot of model H20 slices overlayed with real emissions
plt.scatter(disp_dict['time'][tidx], mod_val)
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 0,0.000002, facecolor='gray', alpha=0.1, edgecolor='w')
plt.ylim([0,0.000002])
plt.xlim([0,2100])
plt.show()

#================================FLIGHT ANIMATION==================================
fig = plt.figure()
ax = p3.Axes3D(fig)

# create initial frame
point, = ax.plot([disp_dict['lcn'][0,0]],[disp_dict['lcn'][0,1]],[disp_dict['lcn'][0,2]], 'o')
ax.contourf(WLAT, WLONG, np.zeros(np.shape(WLAT)), alpha=0.3)
line, = ax.plot(disp_dict['lcn'][:,0], disp_dict['lcn'][:,1], disp_dict['lcn'][:,2], label='flight path', color='gray', alpha=0.3)
ax.legend()
ax.set_xlim([min(disp_dict['lcn'][:,0]), max(disp_dict['lcn'][:,0])])
ax.set_ylim([min(disp_dict['lcn'][:,1]), max(disp_dict['lcn'][:,1])])
ax.set_zlim([min(disp_dict['lcn'][:,2]), max(disp_dict['lcn'][:,2])])
time_text = ax.text(0.05,0.05,0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

#make a list of all times within plume from emissions
smoky = []
for item in emis_dict['smoke']:
	smoky.extend(np.arange(item[0],item[1]))

# move the point position at every frame
def update_point(n, disp_dict,smoky,point):
    point.set_data(np.array([disp_dict['lcn'][n,0],disp_dict['lcn'][n,1]]))
    point.set_3d_properties(disp_dict['lcn'][n,2], 'z')
    time_text.set_text('Time (sec) = %s' %(n*dt))
    if disp_dict['time'][n] in smoky:
    	point.set_color('r')
    else:
    	point.set_color('k')
    return point, time_text,

#plot the first 1500 frames (3000sec) - roughtly the length of the simulation
ani=animation.FuncAnimation(fig, update_point, 1500, fargs=(disp_dict,smoky,point), interval=1)
# ani.save('./test_ani.gif', writer='imagemagick',fps=120)
plt.show()


#================================VIRTICAL PROFILE==================================
#define start and end of the corkscrew in sec from beginning of simulation
csStart = 2400
csEnd = 2700
plt.title('CO2 PROFILE (from Corskrew flight)')
s = np.argmin(abs(disp_dict['time']-csStart))
f = np.argmin(abs(disp_dict['time']-csEnd))
plt.scatter(disp_dict['CO2'][s:f],disp_dict['lcn'][s:f,2] )
plt.ylim([0,2000])
plt.xlabel('CO2 mixing ratio [ppmv]')
plt.ylabel('height [m]')
plt.tight_layout()
plt.savefig(fig_dir + 'CO2ProfileCorkscrew.pdf')
plt.show()

#create a horizontall averaged plume profile
zave = np.nanmean(np.nanmean(np.nanmean(z,0),1),1) 	#average vertical level height)
plume_profile = np.nansum(np.nansum(qvapor,2),2)
plume_profile[plume_profile<0] = np.nan 			#mask negataives
plt.contourf(plume_profile.T,cmap=plt.cm.cubehelix_r)
cbar = plt.colorbar()
cbar.set_label('H2O mixing ratio')
ax = plt.gca()
ax.set_yticks(np.arange(0,np.shape(plume_profile)[1],5))
ax.set_yticklabels(zave[::5])
plt.xlabel('time [min]')
plt.ylabel('height [m]')
plt.title('EVOLUTION OF SMOKE CONCENTRATION COLUMN')
plt.tight_layout()
plt.savefig(fig_dir + 'SmokeColumn.pdf')
plt.show()


# #plot test slice 
# qtest = qvapor[100,3:,:,:]
# qrot = rotate(qtest, mean_wind, axes=(1,2),reshape=True, mode='constant', cval=np.nan, prefilter=True)
# qtot = np.nansum(qrot,2)
# plt.contourf(qtot)
# plt.colorbar()
# plt.show()

# #rotate qvapor grid onto mean wind direction
# rot_QVAPOR = rotate(qvapor, mean_wind, axes=(2,3),reshape=True, mode='constant', cval=np.nan, prefilter=True)


