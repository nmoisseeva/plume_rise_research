# referencing simulation with flight data


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
wrfinterp = '/Users/nadya2/data/plume/RxCADRE/interp/wrfinterp_LG2'
bounds_shape = '/Users/nadya2/data/qgis/LG2012_WGS'
disp_data = '/Users/nadya2/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
emis_data = '/Users/nadya2/data/RxCADRE/dispersion/Data/Emissions_L2G_20121110.csv'

ll_utm = np.array([521620,3376766]) 	#lower left corner of the domain in utm
basemap_path = '/Users/nadya2/code/plume/RxCADRE/npy/%s_%s_bm.npy' %(ll_utm[0],ll_utm[1])
lvl = np.arange(0,2000,50) 				#
emis_excl = 4 							#number of samples to excluded (from the END!)
#=================end of input===============


print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
nc_interp = netcdf.netcdf_file(wrfinterp, mode ='r')  
nc_inputdata = netcdf.netcdf_file(wrfinput, mode ='r')

#get dimensions of the data
nT,nY,nX = np.shape(nc_inputdata.variables['XLONG'])

#create a UTM grid
UTMx = nc_inputdata.variables['XLONG'][0,:,:] + ll_utm[0]
UTMy = nc_inputdata.variables['XLAT'][0,:,:] + ll_utm[1]

#convert coordinate systems to something basemaps can read
wgs84=pyproj.Proj("+init=EPSG:4326")
epsg26916=pyproj.Proj("+init=EPSG:26916")
WGSx, WGSy= pyproj.transform(epsg26916,wgs84,UTMx.ravel(),UTMy.ravel())
WLONG, WLAT = np.reshape(WGSx, np.shape(UTMx)), np.reshape(WGSy, np.shape(UTMy))

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
polygons = bm.readshapefile(bounds_shape,name='fire_bounds',drawbounds=True)
fireim = nc_data.variables['GRNHFX'][25,:,:]
bm.imshow(fireim)
plt.show()

#extract model time info
runstart = nc_data.START_DATE[-8:]
tsec = nc_data.variables['XTIME'][:] * 60 		#get time in seconds since run start
model_ssm = int(runstart[0:2])*3600 + int(runstart[3:5])*60


#================================DISPERSION==================================
#extract and format dispersion data
disp_dict = {}
disp_array = np.genfromtxt(disp_data, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')

start_idx = np.argmin(abs(disp_array[:,0] - model_ssm))
disp_dict['time']= disp_array[start_idx:,0] - model_ssm +1
disp_dict['time'] = disp_dict['time'].astype(int)
disp_dict['CO'] = disp_array[start_idx:,1]
disp_dict['CO2'] = disp_array[start_idx:,2]
disp_dict['CH4'] = disp_array[start_idx:,3]
disp_dict['H2O'] = disp_array[start_idx:,4]
disp_dict['lcn'] = np.array(zip(disp_array[start_idx:,5],disp_array[start_idx:,6],disp_array[start_idx:,7]))
# disp_dict['lcn'] = zip(disp_array[start_idx:,5],disp_array[start_idx:,6])
disp_dict['meta']= 'time: seconds since model start run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation MSL'


#construct KDtree from idealized grid
# Z = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:])/9.81 
tidx = [np.argmin(abs(disp_dict['time']-t)) for t in tsec]
dt = disp_dict['time'][1] - disp_dict['time'][0]

# grid_coord = zip(WGSy,WGSx)
grid_coord = zip(WGSy,WGSx,lvl)
gridTree = KDTree(grid_coord)
dist, grid_id = gridTree.query(np.array(disp_dict['lcn'])[tidx])

#calculate H20 MR departure and extract point locations
qvapor = nc_interp.variables['QVAPOR'][:,:,:,:] - nc_interp.variables['QVAPOR'][0,:,:,:]
mod_val = np.empty((len(tsec))) * np.nan
obs_val = np.empty((len(tsec))) * np.nan
for nt in range(len(tsec)):
	flatq = qvapor[nt,::-1,:,:].ravel()
	mod_val[nt] = flatq[grid_id[nt]]
	mv = disp_dict['H2O'][tidx[nt]]		#in percent by volume
	obs_val[nt] = 0.622*mv/(100-mv)

# plt.plot(mod_val)
# plt.plot(obs_val)
# plt.show()

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



#plot of CO2 slices
plt.figure()
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
plt.scatter(disp_dict['time'],disp_dict['H2O'])
ax = plt.gca()
for nSlice in range(len(emis_dict['smoke'])):
	shade = np.arange(emis_dict['smoke'][nSlice][0],emis_dict['smoke'][nSlice][1])
	ax.fill_between(shade, 0,1.5, facecolor='gray', alpha=0.1, edgecolor='w')
plt.ylim([0,1.5])
plt.xlim([0,3000])
plt.show()


#================================FLIGHT INFO==================================
fig = plt.figure()
ax = p3.Axes3D(fig)

# create the first plot
point, = ax.plot([disp_dict['lcn'][0,0]],[disp_dict['lcn'][0,1]],[disp_dict['lcn'][0,2]], 'o')
ax.plot3D(WLAT[0,:],WLONG[0,:],0, 'k--', lw=2)
ax.plot3D(WLAT[-1,:],WLONG[-1,:],0, 'k--', lw=2)
ax.plot3D(WLAT[:,0],WLONG[:,0],0, 'k--', lw=2)
ax.plot3D(WLAT[:,-1],WLONG[:,-1],0, 'k--', lw=2)

# line, = ax.plot([disp_dict['lcn'][:][0]], [disp_dict['lcn'][:][1]], [disp_dict['lcn'][:][2]], label='parametric curve', color='blue', alpha=0.5)
ax.legend()
ax.set_xlim([min(disp_dict['lcn'][:,0]), max(disp_dict['lcn'][:,0])])
ax.set_ylim([min(disp_dict['lcn'][:,1]), max(disp_dict['lcn'][:,1])])
ax.set_zlim([min(disp_dict['lcn'][:,2]), max(disp_dict['lcn'][:,2])])
time_text = ax.text(0.05,0.05,0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

#make a list of all times within plume from emissions
smoky = []
for item in emis_dict['smoke']:
	smoky.extend(np.arange(item[0],item[1]))


# second option - move the point position at every frame
def update_point(n, disp_dict,smoky,point):
    point.set_data(np.array([disp_dict['lcn'][n,0],disp_dict['lcn'][n,1]]))
    point.set_3d_properties(disp_dict['lcn'][n,2], 'z')
    time_text.set_text('Time (sec) = %s' %(n*dt))
    if disp_dict['time'][n] in smoky:
    	point.set_color('r')
    else:
    	point.set_color('k')
    return point, time_text,

#plot the first 2200 (~35min) - corresponding to length of simulation
ani=animation.FuncAnimation(fig, update_point, 2200, fargs=(disp_dict,smoky,point), interval=1)
# ani.save('./test_ani.gif', writer='imagemagick',fps=120)
plt.show()