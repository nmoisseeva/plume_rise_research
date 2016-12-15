
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

wrfdata = wrfpath +'wrfout_ML'
# wrfinterp = wrfpath + part + '/interp/wrfinterp_' + test
print('Extracting NetCDF data from %s ' %wrfdata)
nc_data = netcdf.netcdf_file(wrfdata, mode ='r')  
xstep = int(nc_data.DX)

#get PT data, mask and save for animation
T = np.copy(nc_data.variables['T'][:,:,:,:])
TKE = np.copy(nc_data.variables['TKE'][:,:,:,:])
tdim,zdim,ydim,xdim = np.shape(T)
slice_y = int(ydim/2)
slice_x = int(xdim/2)

#get heights
print('...calculating height levels')
lvls = (nc_data.variables['PHB'][:,:,:,:] + nc_data.variables['PH'][:,:,:,:])/9.81

nc_data.close()

#==========================plotting===========================

#configure colormap
C = np.genfromtxt('./Tcmap',usecols=(0,1,2))
cm = mpl.colors.ListedColormap(C/255.0)

#animation of TKE growth - manual
fig = plt.figure(figsize=(9,12))
#set up axes
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
ax1.set_title('BL GROWTH', fontsize=13)
ax2.set_title('TKE', fontsize=13)

#animate plots
sub_frm = []
for ti in range(tdim):
	#plot theta
	sub1 = ax1.pcolormesh(T[ti,:,slice_y,:]+300)
	ax1.autoscale_view()
	ax1.set_yticks(np.arange(0,zdim))
	ax1.set_yticklabels(lvls[ti,:,slice_y,slice_x])
	#plot tke
	sub2 = ax2.pcolormesh(TKE[ti,:,slice_y,:], vmin=0, vmax=0.8, cmap=plt.cm.cubehelix_r)
	ax2.set_yticks(np.arange(0,zdim))
	ax2.set_yticklabels(lvls[ti,:,slice_y,slice_x])
	ttl = ax2.annotate('t = %s min' %(ti), xy=(0.45, -0.2), xycoords='axes fraction',fontsize=14, fontweight='bold')
	sub_frm.append([sub1,sub2,ttl])
ani = animation.ArtistAnimation(fig, sub_frm, interval=120, blit=False)
#cosmetics and saving
axes = [ax1, ax2]
# sub1.colorbar()
# sub2.colorbar()
for nAx in axes:
	# nAx.set_yticks(np.arange(0,zdim,10))
	# nAx.set_yticklabels(lvl_hgt[::10])
	nAx.set_ylabel('MSL height [m]')
	nAx.set_xlim([0,xdim])
	nAx.set_ylim([0,zdim])
	nAx.set_xticks(np.arange(0,xdim,20))
	nAx.set_xticklabels(np.arange(0,xdim*xstep,20*xstep))
# fig.subplots_adjust(right=0.85, bottom=0.07, left=0.1,top=0.92)
# fig.colorbar(sub1,cax=cax,label='TKE [m2/s2]')
fig_path = fig_dir +'/BLgrowth.gif'
# ani.save(fig_path, writer='imagemagick',fps=120)
print('.....TKE animation saved as: %s' %fig_path)
plt.show()
# plt.close()

