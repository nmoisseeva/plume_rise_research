
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import matplotlib.animation as animation
import numpy.ma as ma
import matplotlib as mpl 



wrfpath = '/Users/nadya2/data/plume/comps/'
fig_dir = '/Users/nadya2/data/plume/figs/comps/'
part = 'A'
section = ['a_Aug2', 'a_Aug2', 'a_Aug2']
sec_tags = ['CALM', 'LIGHT WINDS', 'MODERATE WINDS']
rx1 = [30,60] 	#averaging grids for ridge 1
rx2 = [90,120]	#averaging grids for ridge 2
fx = [5,20] 		#averaging grids for "flat area"
t_ave = 30 			#averaging window for BL top
profile_data = {'r1':{'tag':'Ridge 1', 'xi':30, 'xf':60}, \
				'r2':{'tag':'Ridge 2', 'xi':90, 'xf':120},\
				'f':{'tag':'flat', 'xi':5, 'xf':20}}


#-----------------------end of input-------------------------

plot_data = {}
for nTest,test in enumerate(section):
	plot_data[test] = {}
	plot_data[test]['tag'] = sec_tags[nTest]
	wrfdata = wrfpath + part + '/wrfout_' + test
	wrfinterp = wrfpath + part + '/interp/wrfinterp_' + test
	print('Extracting NetCDF data from %s ' %wrfdata)
	nc_data = netcdf.netcdf_file(wrfdata, mode ='r')    

	# #sanity check: level height data (otherwise interpolated with ncl)
	# print('Calculating vertical level heights for the domain')
	# lvlhgt = (nc_data.variables['PH'][:,:,:,:] + nc_data.variables['PHB'][:,:,:,:]) / 9.81
	# nc_data.variables['LVLHGT'] = lvlhgt
	# dims = np.shape(lvlhgt)

	#get PT data, mask and save for animation
	interp_data = netcdf.netcdf_file(wrfinterp, mode ='r') 
	interpT = np.copy(interp_data.variables['T'][:,:,:,:])
	tdim,zdim,ydim,xdim = np.shape(interpT)
	interpT[interpT>100] = np.nan
	interp_masked = ma.masked_invalid(interpT)
	plot_data[test]['animation'] = interp_masked

	#create spacially averaged profiles for specified locations
	plot_data[test]['BLheight'], plot_data[test]['BLtemp'] = {}, {}
	for location in profile_data:
		xi = profile_data[location]['xi']
		xf = profile_data[location]['xf']
		profile = np.nanmean(interpT[:,:,75,xi:xf],2)

 		#find height and strength of bl top
 		delT = np.empty((tdim,zdim-1))
 		BLz = np.empty((tdim))
 		BLdT = np.empty((tdim))
 	
		for nTime in range(tdim):
			delT[nTime,:] = profile[nTime,1:] - profile[nTime,0:-1]
			BLz[nTime] = np.argmax(abs(delT[nTime,50:]))
			BLdT[nTime] = delT[nTime,BLz[nTime]]

 		# do temporal averaging
	 	weights = np.repeat(1.0, t_ave)/t_ave
	 	BLH = np.convolve(BLz,weights)
	 	BLT = np.convolve(BLdT,weights)
	 	plot_data[test]['BLheight'][location] = BLH
	 	plot_data[test]['BLtemp'][location] = BLT

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
	plt.legend()
plt.suptitle('BOUNDARY LAYER HEIGHT')
plt.show()








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
ax1.set_xlim([0,xdim])
ax1.set_ylim([0,zdim])
ax2.set_xlim([0,xdim])
ax2.set_ylim([0,zdim])
ax3.set_xlim([0,xdim])
ax3.set_ylim([0,zdim])
#animate plots
sub_frm = []
for ti in range(tdim):
	sub1 = ax1.pcolormesh(plot_data['a_Aug2']['animation'][ti,:,75,:]+300, vmin=300, vmax=310, cmap=cm)
	ax1.autoscale_view()
	sub2 = ax2.pcolormesh(plot_data['a_Aug2']['animation'][ti,:,75,:]+300, vmin=300, vmax=310, cmap=cm)
	sub3 = ax3.pcolormesh(plot_data['a_Aug2']['animation'][ti,:,75,:]+300, vmin=300, vmax=310, cmap=cm)
	sub_frm.append([sub1,sub2,sub3])
ani = animation.ArtistAnimation(fig, sub_frm, interval=120, blit=False)
#cosmetics and saving
fig.subplots_adjust(right=0.85, bottom=0.05, left=0.08,top=0.92)
cax = fig.add_axes([0.88, 0.1, 0.02, 0.8]) #[left, bottom, width, height]
fig.colorbar(sub1,cax=cax,label='Potential temperature [K]')
fig_path = fig_dir + part + '/BL_animation.gif'
ani.save(fig_path)
plt.close()
print('.....BL growth animation saved as: %s' %fig_path)



