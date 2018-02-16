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
# tag = ['W4S400F3R0','W4S400F13R0']
tag = ['W4S400F3R0']

fig_dir = '/Users/nmoisseeva/code/plume/main/figs/'
lvl = np.arange(0,2500,40)	 			#vertical levels in m
# half_width = 20 						#number of grids to use for averaging ALONG fireline
# window = [10,65] 						#averaging moving window - number of grids before and after peak flux

#=================end of input===============

print('ANALYSIS OF PLUME TILT')
print('===================================')


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
	wave_plume = avedict['W'].copy()
	wave_plume[avedict['QVAPOR']<0.1] = np.nan 		#mask where there is no plume
	wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
	wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
	qmax_profile = np.nanmax(avedict['QVAPOR'],1) 	#get max q profile
	qmax_idx = np.nanargmax(avedict['QVAPOR'][np.isfinite(qmax_profile)],1)		#get donwind location
	uave_plume = avedict['U'].copy()
	uave_plume[avedict['QVAPOR']<0.1] = np.nan 		#mask where there is no plume
	umax_profile = np.nanmax(uave_plume,1) 		#get the profiles
	umax_idx = np.nanargmax(uave_plume[np.isfinite(umax_profile)],1)		#get downwind location (index)
	umin_profile = np.nanmin(uave_plume,1) 		#get the profiles
	umin_idx = np.nanargmin(uave_plume[np.isfinite(umax_profile)],1)		#get downwind location (index)
	# uatw_profile = [avedict['U'][ni,i] for ni, i in enumerate(wmax_idx)]	#get u along the profile of maximum w (out of curiosity only)

	#create a vertical slize based at a single downwind location of maximum vertical lift and max q
	qslicew_idx = np.max(wmax_idx)
	qslicew = avedict['QVAPOR'][:,qslicew_idx]
	qsliceq_idx = np.max(qmax_idx)
	qsliceq = avedict['QVAPOR'][:,qsliceq_idx]

	#===========================plotting===========================
	#vertical concentration slice at donwind locations of wmax and qmax
	plt.figure(figsize=(18,6))
	plt.subplot(1,3,1)
	plt.title('VALUE ALONG MAX PROFILE OF W, MAX AND MIN U')
	plt.plot(wmax_profile,lvl,'.-',label='$w_{max}$')
	plt.plot(umax_profile,lvl,'.-',label='$u_{max}$')
	plt.plot(umin_profile,lvl,'.-',label='$u_{min}$')
	# plt.plot(uatw_profile,lvl[np.isfinite(wmax_profile)],'.-',label='$u_{wmax}$')
	plt.xlabel('velocity [m/s]')
	plt.ylabel('height [m]')
	plt.legend()
	plt.ylim([0,lvl[-1]])
	plt.subplot(1,3,2)
	plt.title('HORIZONTAL LOCATION OF EXTREMA')
	plt.plot(wmax_idx,lvl[np.isfinite(wmax_profile)],'.-',label='$w_{max}$')
	plt.plot(qmax_idx,lvl[np.isfinite(qmax_profile)],'.--',label='$q_{max}$')
	plt.plot(umax_idx,lvl[np.isfinite(umax_profile)],'.--',label='$u_{max}$')
	plt.plot(umin_idx,lvl[np.isfinite(umin_profile)],'.--',label='$u_{min}$')
	plt.xlabel('x-grid [#]')
	plt.ylabel('height [m]')
	plt.ylim([0,lvl[-1]])
	plt.legend()
	plt.subplot(1,3,3)
	plt.title('Q CONCENTRATION DOWNWIND (SLICE)')
	plt.plot(qslicew, lvl, '.-',label='based on $w_{max}$')
	plt.plot(qsliceq, lvl, '.--', label='based on $q_{max}$')
	plt.xlabel('water vapor [g/kg]')
	plt.ylabel('height [m]')
	plt.ylim([0,lvl[-1]])
	plt.legend()
	# plt.savefig(fig_dir + 'profiles_%s.pdf' %Case)
	plt.show()
	plt.close()

#
