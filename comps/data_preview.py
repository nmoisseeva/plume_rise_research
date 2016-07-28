
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import matplotlib.animation as animation


wrfdata = '/Users/nadya2/data/plume/comps/part1/wrfout_3'
vert_levels = np.arange(0,1500,50)
wrfinterp = '/Users/nadya2/code/plume/comps/interpT.nc'

# vars_4d = ['U','V','QVAPOR']
# vars_3d = ['HGT','GRNHFX','AVG_FUEL_FRAC','GRNQFX']

print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')             
# fc_data = {}

# #copy variables of interest into a dictionary
# for var4 in vars_4d:
# 	print('......... %s' %var4)
# 	fc_data[var4] = nc_data.variables[var4][:,:,:,:]

# for var3 in vars_3d:
#     print('......... %s' %var3)
#     fc_data[var3] = nc_data.variables[var3][:,:,:]

# # ===========vertical interpolation done in ncl===================
#add level height data
print('Calculating vertical level heights for the domain')
# nc_data.createVariable('LVLHGT', float, ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
nc_data.variables['LVLHGT'] = lvlhgt
dims = np.shape(lvlhgt)


# # #interpolate to vertical height levels
# # print('Performing vertical interpolation of model data')
# # T_agl = np.empty((dims[0],len(vert_levels),dims[2],dims[3]))
# # for nTime in range(dims[0]):
# # 	print('.....timestep since simulation start: %smin' %nc_data.variables['XTIME'][nTime])
# # 	for nY in range(dims[2]):
# # 		for nX in range(dims[3]):
# # 			eta =  lvlhgt[nTime,:-1,nY,nX]
# # 			interp_func = interpolate.interp1d(eta, nc_data.variables['T'][nTime,:,nY,nX],\
# # 						 'cubic', bounds_error=False, fill_value=np.nan)
# # 			T_agl[nTime,:,nY,nX] = interp_func(vert_levels)
# #===========vertical interpolation done in ncl===================

ttest= 170

print("plotting....")
plt.figure()
plt.pcolormesh(nc_data.variables['T'][ttest,0:20,75,:]+300,vmin=300, vmax=310)
plt.colorbar()
plt.show()
plt.close()

plt.figure()
nc_interp_data = netcdf.netcdf_file(wrfinterp, mode ='r') 
interpT = np.copy(nc_interp_data.variables['T'][:,:,:,:])
interpT[interpT>100] = np.nan
plt.pcolormesh(interpT[ttest,:,75,:]+300,vmin=300, vmax=310)
plt.colorbar()
plt.show()
plt.close()


#sanity check: plot vertical column of temperatures for raw and interp dataset
plt.figure()
plt.subplot(2,1,1)
plt.plot(lvlhgt[ttest,:23,75,75],nc_data.variables['T'][ttest,:23,75,75]+300,'k')
plt.subplot(2,1,2)
plt.plot(np.arange(0,1.5,0.02)[3:]*1000,interpT[ttest,:,75,75][3:]+300,'r')
plt.show()
plt.close()


#animation of BL growth
fig = plt.figure()
ims = []
for i in range(len(interpT[:,0,75,0])):
    t_step = int(i)
    im = plt.pcolormesh(interpT[i,:,75,:]+300, vmin=300, vmax=310)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=120, blit=False,repeat_delay=1000)
plt.colorbar()

plt.show()
