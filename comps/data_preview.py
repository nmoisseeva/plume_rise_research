
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



wrfpath = '/Users/nadya2/data/plume/comps/'
fig_dir = '/Users/nadya2/code/plume/figs/comps/'
part = 'D'
section = ['a', 'b', 'c']
sec_tags = ['CALM', 'LIGHT WINDS', 'MODERATE WINDS']
t_ave = 45 										#averaging window for BL top
interp_lvl = [0,2000] 						#interpolation intervals from ncl (m) [min,max]
interp_step = 30
BLi = 500 										#intial residual level height (m)
profile_data = {'r1':{'tag':'Ridge 1', 'xi':55, 'xf':70}, \
				'r2':{'tag':'Ridge 2', 'xi':95, 'xf':110},\
				'f':{'tag':'flat', 'xi':5, 'xf':20}}


#-----------------------end of input-------------------------

print('========COMPS DATA ANALYSIS==========')
print('Performing anlsysis for: PART %s' %part)
print('=====================================')

plot_data = {}

#get initial residual level
lvl_hgt = np.arange(interp_lvl[0],interp_lvl[1],interp_step)
zi = np.argmin(abs(lvl_hgt - BLi))

for nTest,test in enumerate(section):
	plot_data[test] = {}
	plot_data[test]['tag'] = sec_tags[nTest]
	wrfdata = wrfpath + part + '/wrfout_' + test
	wrfinterp = wrfpath + part + '/interp/wrfinterp_' + test
	print('Extracting NetCDF data from %s ' %wrfdata)
	nc_data = netcdf.netcdf_file(wrfdata, mode ='r')   

	#get PT data, mask and save for animation
	interp_data = netcdf.netcdf_file(wrfinterp, mode ='r') 
	interpT = np.copy(interp_data.variables['T'][:,:,:,:])
	interpTKE = np.copy(interp_data.variables['TKE'][:,:,:,:])
	tdim,zdim,ydim,xdim = np.shape(interpT)
	slice_y = int(ydim/2)
	interpT[interpT>100], interpTKE[interpTKE>100] = np.nan, np.nan
	interpT_masked, interpTKE_masked = ma.masked_invalid(interpT), ma.masked_invalid(interpTKE)
	plot_data[test]['Tanim'] = interpT_masked
	plot_data[test]['TKEanim'] = interpTKE_masked
	terrain = nc_data.variables['HGT'][0,slice_y,:]
	terrain = terrain/interp_step
	print('Calculating cross-wind averages')
	plot_data[test]['cwTKE'] = np.nanmean(plot_data[test]['TKEanim'][:,:,:,:],2)
	plot_data[test]['cwT'] = np.nanmean(plot_data[test]['Tanim'][:,:,:,:],2)

	#create spacially averaged profiles for specified locations
	plot_data[test]['BLheight'], plot_data[test]['BLtemp'] = {}, {}
	plot_data[test]['cwBLheight'], plot_data[test]['cwBLtemp'] = {}, {}
	for location in profile_data:
		xi = profile_data[location]['xi']
		xf = profile_data[location]['xf']
		profile = np.nanmean(interpT[:,:,slice_y,xi:xf],2)
		cw_profile = np.nanmean(np.nanmean(interpT[:,:,:,xi:xf],2),2)

 		#find height and strength of bl top
 		delT = np.empty((tdim,zdim-1))
 		BLz = np.empty((tdim))
 		BLdT = np.empty((tdim))
 		cw_delT = np.empty((tdim,zdim-1))
 		cwBLz = np.empty((tdim))
 		cwBLdT = np.empty((tdim))
 	
		for nTime in range(tdim):
			delT[nTime,:] = profile[nTime,1:] - profile[nTime,0:-1]
			BLz[nTime] = np.argmax(abs(delT[nTime,zi:]))
			BLdT[nTime] = delT[nTime,zi+BLz[nTime]]
			cw_delT[nTime,:] = cw_profile[nTime,1:] - cw_profile[nTime,0:-1]
			cwBLz[nTime] = np.argmax(abs(cw_delT[nTime,zi:]))
			cwBLdT[nTime] = cw_delT[nTime,zi+cwBLz[nTime]]

 		# do temporal averaging
	 	weights = np.hanning(t_ave)/sum(np.hanning(t_ave)) 	#create a tapered averaging window
	 	BLH = np.convolve(BLz,weights)
	 	BLT = np.convolve(BLdT,weights)
	 	cwBLH = np.convolve(cwBLz,weights)
	 	cwBLT = np.convolve(cwBLdT,weights)
	 	plot_data[test]['BLheight'][location] = BLH
	 	plot_data[test]['BLtemp'][location] = BLT
	 	plot_data[test]['cwBLheight'][location] = cwBLH
	 	plot_data[test]['cwBLtemp'][location] = cwBLT

	nc_data.close()
	interp_data.close()



#==========================plotting===========================

print('Plotting data:')

fig = plt.figure(figsize=(9,12))
for nTest, test in enumerate(section):
	fig.add_subplot(3,1,nTest+1)
	plt.title(plot_data[test]['tag'])
	for loc in plot_data[test]['BLtemp']:
		plt.plot(plot_data[test]['BLheight'][loc],label=loc)
		plt.xlim([0,tdim])
		plt.xlabel('time [min]')
		plt.ylabel('height MSL [m]')
		ax = plt.gca()
		ax.set_yticks(np.arange(0,zdim-zi,10))
		ax.set_yticklabels(np.arange(BLi,interp_lvl[1],10*interp_step))
	plt.legend(loc='lower right')
plt.suptitle('BOUNDARY LAYER HEIGHT| t_ave = %s' %t_ave,fontweight='bold',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.1)
fig_path = fig_dir + part + '/BL_height.pdf'
plt.savefig(fig_path)
# plt.show()
plt.close()
print('.....BL height timeseries saved as: %s' %fig_path)

fig = plt.figure(figsize=(9,12))
for nTest, test in enumerate(section):
	fig.add_subplot(3,1,nTest+1)
	plt.title(plot_data[test]['tag'])
	for loc in plot_data[test]['BLtemp']:
		plt.plot(plot_data[test]['BLtemp'][loc],label=loc)
		plt.xlim([0,tdim])
		plt.xlabel('time [min]')
		plt.ylabel('inversion strength [K/100m]')
		ax = plt.gca()
		ax.set_yticks(np.arange(0,0.5,0.1))
	plt.legend(loc='lower right')
plt.suptitle('BOUNDARY LAYER STRENGTH | t_ave = %s' %t_ave,fontweight='bold',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.1)
fig_path = fig_dir + part + '/BL_strength.pdf'
plt.savefig(fig_path)
# plt.show()
plt.close()
print('.....BL strength timeseries saved as: %s' %fig_path)



fig = plt.figure(figsize=(9,12))
for nTest, test in enumerate(section):
	fig.add_subplot(3,1,nTest+1)
	plt.title(plot_data[test]['tag'])
	for loc in plot_data[test]['cwBLtemp']:
		plt.plot(plot_data[test]['cwBLheight'][loc],label=loc)
		plt.xlim([0,tdim])
		plt.xlabel('time [min]')
		plt.ylabel('height MSL [m]')
		ax = plt.gca()
		ax.set_yticks(np.arange(0,zdim-zi,10))
		ax.set_yticklabels(np.arange(BLi,interp_lvl[1],10*interp_step))
	plt.legend(loc='lower right')
plt.suptitle('CROSSWIND BOUNDARY LAYER HEIGHT| t_ave = %s' %t_ave,fontweight='bold',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.1)
fig_path = fig_dir + part + '/cwBL_height.pdf'
plt.savefig(fig_path)
# plt.show()
plt.close()
print('.....Crosswind BL height timeseries saved as: %s' %fig_path)

fig = plt.figure(figsize=(9,12))
for nTest, test in enumerate(section):
	fig.add_subplot(3,1,nTest+1)
	plt.title(plot_data[test]['tag'])
	for loc in plot_data[test]['cwBLtemp']:
		plt.plot(plot_data[test]['cwBLtemp'][loc],label=loc)
		plt.xlim([0,tdim])
		plt.xlabel('time [min]')
		plt.ylabel('inversion strength [K/100m]')
		ax = plt.gca()
		ax.set_yticks(np.arange(0,0.5,0.1))
	plt.legend(loc='lower right')
plt.suptitle('CROSSWIND BOUNDARY LAYER STRENGTH | t_ave = %s' %t_ave,fontweight='bold',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.1)
fig_path = fig_dir + part + '/cwBL_strength.pdf'
plt.savefig(fig_path)
# plt.show()
plt.close()
print('.....Crosswind BL strength timeseries saved as: %s' %fig_path)





#time averaged TKE
fig = plt.figure(figsize=(9,12))
for nTest, test in enumerate(section):
	fig.add_subplot(3,1,nTest+1)
	plt.title(plot_data[test]['tag'])
	# plt.pcolormesh(np.nanmean(plot_data[test]['TKEanim'][:,:,slice_y,:],0), vmin=0, vmax=1, cmap=plt.cm.cubehelix_r)
	pl = plt.pcolormesh(np.nanmean(plot_data[test]['cwTKE'],0), vmin=0, vmax=0.8, cmap=plt.cm.cubehelix_r)
	# plt.plot(np.arange(0,xdim),terrain,c='w',linewidth=20)
	ax = plt.gca()
	ax.set_yticks(np.arange(0,zdim,10))
	ax.set_yticklabels(lvl_hgt[::10])
	ax.set_ylabel('MSL height [m]')
	ax.set_xlim([0,xdim])
	ax.set_ylim([0,zdim])
plt.suptitle('CROSSWIND TIME-AVERAGED TKE',fontweight='bold',fontsize=16)
plt.tight_layout()
fig.subplots_adjust(right=0.85, bottom=0.07, left=0.1,top=0.92)
cax = fig.add_axes([0.88, 0.1, 0.02, 0.8]) #[left, bottom, width, height]
fig.colorbar(pl,cax=cax,label='TKE [m2/s2]')
fig_path = fig_dir + part + '/TKE_cw_ave.pdf'
plt.savefig(fig_path)
# plt.show()
plt.close()
print('.....Crosswind time-averaged TKE saved as: %s' %fig_path)



#animation of BL growth - manual
fig = plt.figure(figsize=(9,12))
plt.suptitle('BL GROWTH | PART %s ' %part, fontsize=16, fontweight='bold')
#configure colormap
C = np.genfromtxt('./Tcmap',usecols=(0,1,2))
cm = mpl.colors.ListedColormap(C/255.0)
#set up axes
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)
ax1.set_title('CALM', fontsize=13)
ax2.set_title('LIGHT WINDS (3 m/s)', fontsize=13)
ax3.set_title('MODERATE WINDS (6 m/s)', fontsize=13)
#animate plots
sub_frm = []
for ti in range(tdim):
	sub1 = ax1.pcolormesh(plot_data[section[0]]['Tanim'][ti,:,slice_y,:]+300, vmin=300, vmax=310, cmap=cm)
	ax1.autoscale_view()
	sub2 = ax2.pcolormesh(plot_data[section[1]]['Tanim'][ti,:,slice_y,:]+300, vmin=300, vmax=310, cmap=cm)
	sub3 = ax3.pcolormesh(plot_data[section[2]]['Tanim'][ti,:,slice_y,:]+300, vmin=300, vmax=310, cmap=cm)
	ttl = ax3.annotate('t = %s min' %(ti), xy=(0.5, -0.2), xycoords='axes fraction',fontsize=14, fontweight='bold')
	sub_frm.append([sub1,sub2,sub3,ttl])
ani = animation.ArtistAnimation(fig, sub_frm, interval=120, blit=False)
#cosmetics and saving
axes = [ax1, ax2, ax3]
for nAx in axes:
	# nAx.plot(np.arange(0,xdim),terrain,c='w',linewidth=10)
	nAx.set_yticks(np.arange(0,zdim,10))
	nAx.set_yticklabels(lvl_hgt[::10])
	nAx.set_ylabel('MSL height [m]')
	nAx.set_xlim([0,xdim])
	nAx.set_ylim([0,zdim])
fig.subplots_adjust(right=0.85, bottom=0.07, left=0.1,top=0.92)
cax = fig.add_axes([0.88, 0.1, 0.02, 0.8]) #[left, bottom, width, height]
fig.colorbar(sub1,cax=cax,label='Potential temperature [K]')
fig_path = fig_dir + part + '/BL_animation.gif'
ani.save(fig_path, writer='imagemagick',fps=120)
print('.....BL growth animation saved as: %s' %fig_path)
plt.close()

#animation of TKE growth - manual
fig = plt.figure(figsize=(9,12))
plt.suptitle('TKE | PART %s ' %part, fontsize=16, fontweight='bold')
#set up axes
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)
ax1.set_title('CALM', fontsize=13)
ax2.set_title('LIGHT WINDS (3 m/s)', fontsize=13)
ax3.set_title('MODERATE WINDS (6 m/s)', fontsize=13)
#animate plots
sub_frm = []
for ti in range(tdim):
	sub1 = ax1.pcolormesh(plot_data[section[0]]['TKEanim'][ti,:,slice_y,:], vmin=0, vmax=0.8, cmap=plt.cm.cubehelix_r)
	ax1.autoscale_view()
	sub2 = ax2.pcolormesh(plot_data[section[1]]['TKEanim'][ti,:,slice_y,:], vmin=0, vmax=0.8, cmap=plt.cm.cubehelix_r)
	sub3 = ax3.pcolormesh(plot_data[section[2]]['TKEanim'][ti,:,slice_y,:], vmin=0, vmax=0.8, cmap=plt.cm.cubehelix_r)
	ttl = ax3.annotate('t = %s min' %(ti), xy=(0.5, -0.2), xycoords='axes fraction',fontsize=14, fontweight='bold')
	sub_frm.append([sub1,sub2,sub3,ttl])
ani = animation.ArtistAnimation(fig, sub_frm, interval=120, blit=False)
#cosmetics and saving
axes = [ax1, ax2, ax3]
for nAx in axes:
	# nAx.plot(np.arange(0,xdim),terrain,c='w',linewidth=10)
	nAx.set_yticks(np.arange(0,zdim,10))
	nAx.set_yticklabels(lvl_hgt[::10])
	nAx.set_ylabel('MSL height [m]')
	nAx.set_xlim([0,xdim])
	nAx.set_ylim([0,zdim])
fig.subplots_adjust(right=0.85, bottom=0.07, left=0.1,top=0.92)
cax = fig.add_axes([0.88, 0.1, 0.02, 0.8]) #[left, bottom, width, height]
fig.colorbar(sub1,cax=cax,label='TKE [m2/s2]')
fig_path = fig_dir + part + '/TKE_animation.gif'
ani.save(fig_path, writer='imagemagick',fps=120)
print('.....TKE animation saved as: %s' %fig_path)
plt.close()

