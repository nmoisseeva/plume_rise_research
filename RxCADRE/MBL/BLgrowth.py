
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import matplotlib.animation as animation
import numpy.ma as ma
import matplotlib as mpl 
import warnings
warnings.filterwarnings("ignore")



wrfpath = '/Users/nadya2/data/plume/RxCADRE/MBL/'
fig_dir = '/Users/nadya2/code/plume/figs/RxCADRE/'
t_ave = 45 										#averaging window for BL top
# interp_lvl = [0,2000] 						#interpolation intervals from ncl (m) [min,max]
# interp_step = 30
# BLi = 500 										#intial residual level height (m)


#-----------------------end of input-------------------------

plot_data = {}

wrfdata = wrfpath +'wrfout_MBL'
# wrfinterp = wrfpath + part + '/interp/wrfinterp_' + test
print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
xstep = int(nc_data.DX)


#get PT data
T = np.copy(nc_data.variables['T'][:,:,:,:])
times,zdim = np.shape(T)[0:2]
profile = np.empty((times,zdim)) * np.nan
for nT in range(times):
	pr = T[nT,:,:,:]
	profile[nT,:] = np.mean(np.mean(pr,1),1)+300
	# pr = T[nT,:,125,75]
	# profile[nT,:] = pr+300


nc_data.close()

#==========================plotting===========================


fig = plt.figure(figsize=(9,12))
plt.imshow(profile.T, origin='lower')
plt.colorbar()
plt.show()


#save path 
# fig_path = fig_dir +'/BLgrowth_spinup.pdf'
# print('.....BL growth plot saved as: %s' %fig_path)
# # plt.close()




