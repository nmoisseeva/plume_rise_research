
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import matplotlib.animation as animation


wrfpath = '/Users/nadya2/data/plume/comps/'
fig_dir = 
part = 'A'
section = ['a_Aug2', 'a_Aug2', 'a_Aug2']

#-----------------------end of input-------------------------

plot_data = {'animation':{}}
for nTest,test in enumerate(section):
	wrfdata = wrfpath + part + '/wrfout_' + test
	wrfinterp = wrfpath + part + '/interp/wrfinterp_' + test
	print('Extracting NetCDF data from %s ' %wrfdata)
	nc_data = netcdf.netcdf_file(wrfdata, mode ='r')    

	# #sanity check: level height data (otherwise interpolated with ncl)
	# print('Calculating vertical level heights for the domain')
	# lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
	# nc_data.variables['LVLHGT'] = lvlhgt
	# dims = np.shape(lvlhgt)

	interp_data = netcdf.netcdf_file(wrfinterp, mode ='r') 
	interpT = np.copy(interp_data.variables['T'][:,:,:,:])
	interpT[interpT>100] = np.nan

	plot_data['animation'][test] = interpT
	nc_data.close()
	interp_data.close()



#==========================plotting===========================

#create subplot of animations



# ttest= 170

# print("plotting....")
# plt.figure()
# plt.pcolormesh(nc_data.variables['T'][ttest,0:20,75,:]+300,vmin=300, vmax=310)
# plt.colorbar()
# plt.show()
# plt.close()

# plt.figure()
# nc_interp_data = netcdf.netcdf_file(wrfinterp, mode ='r') 
# interpT = np.copy(nc_interp_data.variables['T'][:,:,:,:])
# interpT[interpT>100] = np.nan
# plt.pcolormesh(interpT[ttest,:,75,:]+300,vmin=300, vmax=310)
# plt.colorbar()
# plt.show()
# plt.close()

# #sanity check: plot vertical column of temperatures for raw and interp dataset
# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(lvlhgt[ttest,:23,75,75],nc_data.variables['T'][ttest,:23,75,75]+300,'k')
# plt.subplot(2,1,2)
# plt.plot(np.arange(0,1.5,0.02)[3:]*1000,interpT[ttest,:,75,75][3:]+300,'r')
# plt.show()
# plt.close()

tdim,zdim,ydim,xdim = np.shape(interpT)

#animation of BL growth - manual
fig = plt.figure(figsize=(9,12))
plt.suptitle('BL GROWTH | PART %s ' %part, fontsize=16, fontweight='bold')
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)

ax1.set_title('CALM', fontsize=13)
ax2.set_title('LIGHT WINDS (3 m/s)', fontsize=13)
ax3.set_title('MODERATE WINDS (6 m/s)', fontsize=13)

ax1.set_xlim([0,xdim])
ax1.set_ylim([0,zdim])
ax2.set_xlim([0,xdim])
ax2.set_ylim([0,zdim])
ax3.set_xlim([0,xdim])
ax3.set_ylim([0,zdim])

sub_frm = []
for ti in range(tdim):
	sub1 = ax1.pcolormesh(plot_data['animation']['a_Aug2'][ti,:,75,:]+300, vmin=300, vmax=310)
	ax1.autoscale_view()
	sub2 = ax2.pcolormesh(plot_data['animation']['a_Aug2'][ti,:,75,:]+300, vmin=300, vmax=310)
	sub3 = ax3.pcolormesh(plot_data['animation']['a_Aug2'][ti,:,75,:]+300, vmin=300, vmax=310)
	sub_frm.append([sub1,sub2,sub3])
ani = animation.ArtistAnimation(fig, sub_frm, interval=120, blit=False)

fig.subplots_adjust(right=0.85, bottom=0.05, left=0.08,top=0.92)
cax = fig.add_axes([0.88, 0.1, 0.02, 0.8]) #[left, bottom, width, height]
fig.colorbar(sub1,cax=cax,label='Potential temperature [K]')

plt.show()
fig_path = fig_dir + part + '/BL_animation.gif'
ani.save(fig_path)






# plt.savefig(fig_path)
# plt.close()
# print('.....QVAPOR departure (individual subplots) saved as: %s' %fig_path)
# nc_data.close()
