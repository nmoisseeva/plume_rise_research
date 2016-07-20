#==============================================================
# Script to produce plots for sensitivity analysis of plume rise
# Creates the following figures:
#			- subplots of individual plumes (horizontally integrated)
#			- scatter plot of sensitivity
#===============================================================


sens_var = 'wind'
tag = 'v'
test_range = range(2,14,2)
data_path = '/Users/nadya2/data/plume/sensitivity/'
fig_dir = '/Users/nadya2/code/plume/figs/sensitivity/'
plot_dims = [2,3]

#----------------------------end of input------------------------

from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf


print('Performing analysis of individual runs for: %s' %sens_var)
row_int_q = []
fig, axes = plt.subplots(nrows=plot_dims[0], ncols=plot_dims[1],figsize=(20,10))
for nTest, test in enumerate(test_range):
	wrfdata = data_path + sens_var + '/wrfout_%s%s' %(test,tag)

	print('.....Extracting NetCDF data from %s ' %wrfdata)
	nc_data = netcdf.netcdf_file(wrfdata, 'r')             

	#create height vector from geopotential
	lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
	height_vector = np.mean(np.mean(np.mean(lvlhgt,3),2),0)

	#calculate vapor mixing ratio anomaly
	anom = np.nansum(np.nansum(nc_data.variables['QVAPOR'][:,:,:,:] - nc_data.variables['QVAPOR'][0,:,:,:],3),2)
	row_int_q.append(anom)
	plt.subplot(2,3,nTest+1)
	im = plt.pcolormesh(row_int_q[nTest].T, vmin=-0.1, vmax=0.3)
	plt.title('%s: %s%s' %(sens_var,test,tag), fontsize = 12)
	plt.xlabel('time since start [min]')
	plt.ylabel('height [m]')
	plt.ylim([0,nc_data.dimensions['bottom_top']])
	plt.yticks(np.arange(0,nc_data.dimensions['bottom_top'],10), height_vector[0::10])
plt.suptitle('SENSTIVITY ANALYSIS: HORIZONTALLY INTEGRATED QVAPOR ANOMALY')
plt.tight_layout()
fig.subplots_adjust(top=0.92, bottom=0.15)
cax = fig.add_axes([0.1, 0.05, 0.8, 0.03])
fig.colorbar(im, cax=cax,label='H20 mixing ratio departure [kg/kg]', orientation='horizontal')
fig_path = fig_dir + sens_var + '/QVAPOR_departure.pdf'
plt.savefig(fig_path)
plt.close()
print('.....QVAPOR departure (individual subplots) saved as: %s' %fig_path)
nc_data.close()

print('Determining plume rise sensitivity')
qvap_array = np.array(row_int_q)
ave_idx = np.argmax(np.mean(qvap_array,1),1)
end_idx = np.argmax(qvap_array[:,-1,:],1)

max_pr_ave = height_vector[ave_idx]
max_pr_end = height_vector[end_idx]

plt.scatter(test_range,max_pr_ave,color='r',label='time average')
plt.scatter(test_range,max_pr_end,color='b',label='end rise')
plt.legend()
plt.title('SENSITIVITY PLOT: %s' %sens_var)
plt.xlabel('wind velocity [m/s]')
plt.ylabel('plume rise [m]')
fig_path = fig_dir + sens_var + '/QVAPOR_scatter.pdf'
plt.savefig(fig_path)
plt.close()
print('.....QVAPOR sensitivity scatter saved as: %s' %fig_path)
