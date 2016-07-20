
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

wrfdata = '/Users/nadya2/data/plume/comps/part1/wrfout_d01_0000-01-01_00:00:00'
vars_4d = ['U','V','QVAPOR']
vars_3d = ['HGT','GRNHFX','AVG_FUEL_FRAC','GRNQFX','PBLH']

print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='a')             
fc_data = {}

#copy variables of interest into a dictionary
for var4 in vars_4d:
	print('......... %s' %var4)
	fc_data[var4] = nc_data.variables[var4][:,:,:,:]

for var3 in vars_3d:
    print('......... %s' %var3)
    fc_data[var3] = nc_data.variables[var3][:,:,:]

print('Appending netcdf file with level height data')
#add level height data
nc_data.createVariable('LVLHGT', float, ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
nc_data.variables['LVLHGT'] = lvlhgt
nc_data.sync()

nc_data.close()
# plt.imshow(fc_data['V'][100,:,75,:], origin='lower')
