
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate

wrfdata = '/Users/nadya2/data/plume/comps/part1/wrfout_2'
vert_levels = np.arange(0,3000,100)

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


#add level height data
print('Calculating vertical level heights for the domain')
nc_data.createVariable('LVLHGT', float, ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
nc_data.variables['LVLHGT'] = lvlhgt
dims = np.shape(lvlhgt)


#interpolate to vertical height levels
print('Performing vertical interpolation of model data')
T_agl = np.empty((dims[0],len(vert_levels),dims[2],dims[3]))
for nTime in range(dims[0]):
	print('.....timestep since simulation start: %smin' %nc_data.variables['XTIME'][nTime])
	for nY in range(dims[2]):
		for nX in range(dims[3]):
			eta =  lvlhgt[nTime,:-1,nY,nX]
			interp_func = interpolate.interp1d(eta, nc_data.variables['T'][nTime,:,nY,nX],\
						 'cubic', bounds_error=False, fill_value=np.nan)
			T_agl[nTime,:,nY,nX] = interp_func(vert_levels)

# nc_data.close()
plt.pcolormesh(T_agl[120,:,75,:])
plt.show()