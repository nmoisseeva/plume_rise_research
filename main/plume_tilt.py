# nmoisseeva@eoas.ubc.ca
# January 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf
import cmocean
import sys

#====================INPUT===================
wrfdir = '/Users/nmoisseeva/data/plume/main/'
#loop through all
tag = ['W4S400F1R0','W4S400F2R0','W4S400F3R0','W4S400F5R0',
		'W4S400F6R0','W4S400F7R0','W4S400F8R0','W4S400F9R0',
		'W4S400F10R0','W4S400F11R0','W4S400F12R0','W4S400F13R0',
		'W2S400F3R0','W6S400F3R0','W8S400F3R0','W10S400F3R0',
		'W4S0F3R0','W4S700F3R0','W4Sn400F3R0']

#loop through fuel
# tag = ['W4S400F1R0','W4S400F2R0','W4S400F3R0','W4S400F6R0','W4S400F7R0',
# 		'W4S400F8R0','W4S400F9R0','W4S400F11R0','W4S400F12R0']
#loop through wind
# tag = ['W2S400F3R0','W4S400F3R0','W6S400F3R0','W8S400F3R0','W10S400F3R0']
# tag = ['W4S400F3R0']

fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'
lvl = np.arange(0,2500,40)	 			#vertical levels in m
dx = 40.
# half_width = 20 						#number of grids to use for averaging ALONG fireline
# window = [10,65] 						#averaging moving window - number of grids before and after peak flux

#=================end of input===============
#tag reading function read_tag(variable type, string array)
def read_tag(str_tag, str_array):
	out_array = []
	for nTag,tag in enumerate(str_array):
		letters = re.split('\d+', tag) 
		numbers = re.findall('\d+', tag)
		if str_tag=='W':
			out_array.append(int(numbers[0]))
		elif str_tag=='S':
			if 'Sn' in letters:
				out_array.append(int(numbers[1])*-1)
			else:
				out_array.append(int(numbers[1]))
		elif str_tag=='F':
			out_array.append(int(numbers[2]))
	out_array = np.array(out_array)
	return out_array











#============================================

print('ANALYSIS OF PLUME TILT')
print('===================================')

param_dict = {'fire':[]}
profile_dict = {'wmax':[],'q_wmax':[],'tilt':[]}
profile_dict['meta'] = 'wmax: profiles of maximum velocity; \
						q_wmax: tracer profile at downwind location of maximum vertical velocity \
						tilt: grid number of wmax profile locations'

# tilt_vars = np.empty((len(tag),4))
plume_tops = np.empty(len(tag))
cumT = np.empty(len(tag))
meanU = np.empty(len(tag))
plumeTilt = np.empty(len(tag))
fireLine = np.empty((len(tag),2)) #intensity, width 

for nCase,Case in enumerate(tag):
	print('Examining case: %s ' %Case)

	#----------check for interpolated data----------------------------
	avepath = wrfdir + 'interp/wrfave_' + Case + '.npy'

	if os.path.isfile(avepath):
		print('Averaged data found at: %s' %avepath)
		avedict = np.load(avepath).item()   # load here the above pickle
	else:
		sys.exit('ERROR: no averaged data found - run prep_plumes.py first!')

	#extract lcoations of max w, q, u, and minimum u
	qmax_profile = np.nanmax(avedict['QVAPOR'],1) 	#get max q profile
	qmax_profile[qmax_profile<0.1] = np.nan
	qmax_idx = np.nanargmax(avedict['QVAPOR'][np.isfinite(qmax_profile)],1)		#get donwind location
	qmax_meters = qmax_idx*dx

	wave_plume = avedict['W'].copy()
	wave_plume[avedict['QVAPOR']<0.1] = np.nan 		#mask where there is no plume
	wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
	wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
	watq_profile = np.array([avedict['W'][ni,i] for ni, i in enumerate(qmax_idx)])	#get the profiles

	#get plume tilt through linear regression of max concentration values
	tilt = np.poly1d(np.polyfit(qmax_idx,lvl[np.isfinite(qmax_profile)],1))
	#get average fireline intensity and wind
	U = np.nanmean(avedict['U'][3:,0]) #exclude bottom three layers due to friction
	ignited = [i for i in avedict['GHFX'] if i > 2]
	totI = sum(np.array(ignited)*dx)


	#store data
	meanU[nCase] = U
	plumeTilt[nCase] = tilt.c[0]
	fireLine[nCase,:] = totI,len(ignited)
	plume_tops[nCase] = lvl[len(wmax_idx)-1]
	print('plume top: %s' %lvl[len(wmax_idx)-1])

	#create a vertical slize based at a single downwind location of maximum plume rise
	qslicewtop = avedict['QVAPOR'][:,wmax_idx[-1]] #This should be the proper definition

	#get "cumulutive" temperature for the profile
	cumT[nCase] = np.sum(avedict['T'][:len(wmax_idx),0]*dx)

	#===========================plotting===========================
	#vertical concentration slice at donwind locations of wmax and qmax
	plt.figure(figsize=(18,6))
	plt.suptitle('%s' %Case)
	plt.subplot(1,3,1)
	plt.title('VALUE ALONG MAX PROFILE OF W, MAX AND MIN U')
	plt.plot(wmax_profile,lvl,'.-',label='$w_{max}$')
	plt.plot(watq_profile,lvl[np.isfinite(qmax_profile)],'k.-',label='$w_{qmax}$')
	plt.xlabel('velocity [m/s]')
	plt.ylabel('height [m]')
	plt.legend()
	plt.ylim([0,lvl[-1]])
	plt.subplot(1,3,2)

	plt.title('HORIZONTAL LOCATION OF EXTREMA')
	plt.plot(wmax_idx,lvl[np.isfinite(wmax_profile)],'.-',label='$w_{max}$')
	plt.plot(qmax_idx,lvl[np.isfinite(qmax_profile)],'k.--',label='$q_{max}$')
	plt.plot(qmax_idx, tilt(qmax_idx))
	plt.xlabel('x distance [m]')
	ax = plt.gca()
	ax.set_xticks(np.arange(0,75,15))
	ax.set_xticklabels(np.arange(0,75,15)*40)
	plt.ylabel('height [m]')
	plt.ylim([0,lvl[-1]])
	plt.legend()
	plt.subplot(1,3,3)

	plt.title('Q CONCENTRATION DOWNWIND (SLICE)')
	plt.plot(qslicewtop, lvl, 'g.--', label='based on $w_{top}$')
	plt.xlabel('water vapor [g/kg]')
	plt.ylabel('height [m]')
	plt.ylim([0,lvl[-1]])
	plt.legend()
	plt.savefig(fig_dir + 'profiles_%s.pdf' %Case)
	# plt.show()
	plt.close()

#-------------subplots of tilt predictors-----------
#create scatter plots of tilts vs other variables
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('FIRELINE INTENSITY vs TILT')
plt.scatter(fireLine[:,0],plumeTilt,c=meanU)
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('plume tilt')

plt.subplot(1,3,2)
plt.title('WIND vs TILT')
plt.scatter(meanU,plumeTilt,c=fireLine[:,0])
plt.xlabel('mean wind (m/s)')
plt.ylabel('plume tilt')

plt.subplot(1,3,3)
plt.title('WIDTH vs TILT')
plt.scatter(fireLine[:,1],plumeTilt,c=meanU)
plt.xlabel('fireline width [#grids]')
plt.ylabel('plume tilt')
plt.show()
plt.close()



plt.title('FIRELINE INTENSITY vs CumT')
ax = plt.gca()
ax.scatter(fireLine[:,0],cumT, c=meanU)
ax.scatter(fireLine[-3:,0],cumT[-3:], c=np.array([0,700,-400]), cmap=plt.cm.PiYG)
# plt.colorbar()
for i, txt in enumerate(read_tag('F',tag)):
    ax.annotate(txt, (fireLine[i,0]+100,cumT[i]+100), fontsize=9)

plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('cumulative temperature [K]')
plt.show()

plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
sc = ax.scatter(fireLine[:,0]/(meanU),cumT, c=read_tag('S',tag), cmap=plt.cm.PiYG)
cbar=plt.colorbar(sc)
cbar.set_label('surface heat flux [$W/m^2$]')
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('cumulative temperature [K]')
plt.show()

