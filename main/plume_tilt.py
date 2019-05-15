# nmoisseeva@eoas.ubc.ca
# January 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
#import wrf
#import cmocean
import sys
import imp

#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

#=================end of input===============

print('ANALYSIS OF PLUME TILT')
print('===================================')


param_dict = {'fire':[]}
profile_dict = {'wmax':[],'q_wmax':[],'tilt':[]}
profile_dict['meta'] = 'wmax: profiles of maximum velocity; \
						q_wmax: tracer profile at downwind location of maximum vertical velocity \
						tilt: grid number of wmax profile locations'

# tilt_vars = np.empty((len(tag),4))
plume_tops = np.empty(len(plume.tag))
cumT = np.empty(len(plume.tag))
meanU = np.empty(len(plume.tag))
plumeTilt = np.empty(len(plume.tag))
fireLine = np.empty((len(plume.tag),2)) #intensity, width
Qprofiles = np.empty((len(plume.tag),len(plume.lvl)))

for nCase,Case in enumerate(plume.tag):
	print('Examining case: %s ' %Case)

	#----------check for interpolated data----------------------------
	avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'

	if os.path.isfile(avepath):
		print('Averaged data found at: %s' %avepath)
		avedict = np.load(avepath).item()   # load here the above pickle
	else:
		sys.exit('ERROR: no averaged data found - run prep_plumes.py first!')

	#extract lcoations of max w, q, u, and minimum u
	qmax_profile = np.nanmax(avedict['qvapor'],1) 	#get max q profile
	qmax_profile[qmax_profile<0.1] = np.nan
	qmax_idx = np.nanargmax(avedict['qvapor'][np.isfinite(qmax_profile)],1)		#get donwind location
	qmax_meters = qmax_idx*plume.dx

	wave_plume = avedict['w'].copy()
	wave_plume[avedict['qvapor']<0.1] = np.nan 		#mask where there is no plume
	wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
	wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
	watq_profile = np.array([avedict['w'][ni,i] for ni, i in enumerate(qmax_idx)])	#get the profiles

	#get plume tilt through linear regression of max concentration values
	tilt = np.poly1d(np.polyfit(qmax_idx,plume.lvl[np.isfinite(qmax_profile)],1))

	#get fireline intensity and wind
	# U = np.nanmean(avedict['u'][3:,0]) #exclude bottom three layers due to friction
	U = avedict['u'][2,0]

	ignited = np.array([i for i in avedict['ghfx'] if i > 2])


	# #calculate total flux from ground and fire
	# fireI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 )
	# grndI = plume.read_tag('S',[Case])[0]* plume.dx * len(ignited)/(1.2 * 1005)
	# totI = fireI + grndI

	totI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 ) #convert to kinematic heat flux (and remove -kilo)

	print('Mean wind speed: %s std %.2f' %(U,np.nanstd(avedict['u'][3:,0])))

	#store data
	meanU[nCase] = U

	plumeTilt[nCase] = tilt.c[0]
	fireLine[nCase,:] = totI,len(ignited)
	plume_tops[nCase] = plume.lvl[len(wmax_idx)-1]
	print('plume top: %s' %plume.lvl[len(wmax_idx)-1])

	#create a vertical slize based at a single downwind location of maximum plume rise
	qslicewtop = avedict['qvapor'][:,wmax_idx[-1]] #This should be the proper definition

	#get "cumulutive" temperature for the profile
	cumT[nCase] = np.sum(avedict['temp'][:len(wmax_idx),0]*plume.dz)

	#get cumulative delT
	# delT = avedict['temp'][1:len(wmax_idx),0]-avedict['temp'][0:len(wmax_idx)-1,0]
	# cumT[nCase] = np.sum(delT * plume.dz)

	#===========================plotting===========================
	#vertical concentration slice at donwind locations of wmax and qmax
	plt.figure(figsize=(18,6))
	plt.suptitle('%s' %Case)
	plt.subplot(1,3,1)
	plt.title('VALUE ALONG MAX PROFILE OF W, MAX AND MIN U')
	plt.plot(wmax_profile,plume.lvl,'.-',label='$w_{max}$')
	plt.plot(watq_profile,plume.lvl[np.isfinite(qmax_profile)],'k.-',label='$w_{qmax}$')
	plt.xlabel('velocity [m/s]')
	plt.ylabel('height [m]')
	plt.legend()
	plt.ylim([0,plume.lvl[-1]])
	plt.subplot(1,3,2)

	plt.title('HORIZONTAL LOCATION OF EXTREMA')
	plt.plot(wmax_idx,plume.lvl[np.isfinite(wmax_profile)],'.-',label='$w_{max}$')
	plt.plot(qmax_idx,plume.lvl[np.isfinite(qmax_profile)],'k.--',label='$q_{max}$')
	plt.plot(qmax_idx, tilt(qmax_idx))
	plt.xlabel('x distance [m]')
	ax = plt.gca()
	ax.set_xticks(np.arange(0,75,15))
	ax.set_xticklabels(np.arange(0,75,15)*40)
	plt.ylabel('height [m]')
	plt.ylim([0,plume.lvl[-1]])
	plt.legend()
	plt.subplot(1,3,3)

	plt.title('Q CONCENTRATION DOWNWIND (SLICE)')
	plt.plot(qslicewtop, plume.lvl, 'g.--', label='based on $w_{top}$')
	plt.xlabel('water vapor [g/kg]')
	plt.ylabel('height [m]')
	plt.ylim([0,plume.lvl[-1]])
	plt.legend()
	plt.savefig(plume.figdir + 'profiles_%s.pdf' %Case)
	# plt.show()
	plt.close()

	Qprofiles[nCase,:] = qslicewtop



#---------------------calculations over all plumes--------------
normI = fireLine[:,0]/(meanU) #normalized fire intensity
print('Variability of normalized intensity:')
print(np.std(normI))

#-------------subplots of tilt predictors-----------
#create scatter plots of tilts vs other variables
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('FIRELINE INTENSITY vs TILT')
ax1 = plt.scatter(fireLine[:,0],plumeTilt,c=meanU,cmap=plt.cm.viridis)
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('plume tilt')
plt.colorbar(ax1,label='windspeed [m/s]')

plt.subplot(1,3,2)
plt.title('WIND vs TILT')
ax2 = plt.scatter(meanU,plumeTilt,c=fireLine[:,0],cmap=plt.cm.plasma)
plt.xlabel('mean wind (m/s)')
plt.ylabel('plume tilt')
plt.colorbar(ax2, label='fireline intensity [kW/m]')

plt.subplot(1,3,3)
plt.title('WIDTH vs TILT')
ax3 = plt.scatter(fireLine[:,1],plumeTilt,c=meanU,cmap=plt.cm.viridis)
plt.xlabel('fireline width [#grids]')
plt.ylabel('plume tilt')
plt.colorbar(ax3,label='windspeed [m/s]')

plt.tight_layout()
plt.savefig(plume.figdir + 'tilt_subplots.pdf')
# plt.show()
plt.close()


plt.title('FIRELINE INTENSITY vs CumT')
ax1 = plt.scatter(fireLine[:,0],cumT, c=meanU)
# ax = plt.gca()
# ax.scatter(fireLine[:,0],cumT, c=meanU)
# ax.scatter(fireLine[-3:,0],cumT[-3:], c=np.array([0,700,-400]), cmap=plt.cm.PiYG)
# for i, txt in enumerate(plume.read_tag('F',plume.tag)):
#     ax.annotate(txt, (fireLine[i,0]+100,cumT[i]+100), fontsize=9)
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('cumulative temperature [K m]')
plt.colorbar(ax1,label='wind speed [m/s]')
plt.tight_layout()
plt.savefig(plume.figdir + 'I_cumT.pdf')
plt.show()
plt.close()



#get linear regression for the normalized fireLine
regF = np.poly1d(np.polyfit(normI,cumT,1))


plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
sc = ax.scatter(normI,cumT, c=plume.read_tag('S',plume.tag), cmap=plt.cm.PiYG_r, vmin=-600, vmax=600)
plt.plot(normI, regF(normI))
# sc = ax.scatter(fireLine[:,0]/(plume.read_tag('W',plume.tag)),cumT, c=plume.read_tag('S',plume.tag), cmap=plt.cm.PiYG)
for i, txt in enumerate(plume.read_tag('W',plume.tag)):
    ax.annotate(txt, (normI[i]+100,cumT[i]+100), fontsize=9)
plt.colorbar(sc, label='surface heat flux [$W/m^{2}$]')
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.tight_layout()
plt.savefig(plume.figdir + 'normI_cumT.pdf')
plt.show()
plt.close()



plt.figure(figsize=(12,6))
clr = plt.cm.plasma(plt.Normalize()(fireLine[:,0]))
# clr[..., -1] = plt.Normalize()(fireLine[:,0])
for nCase,Case in enumerate(plume.tag):
	plt.subplot(1,2,2)
	plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
	plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile ')
	# plt.colorbar(label='fireline intensity ')
	plt.ylim([0,2])
	plt.subplot(1,2,1)
	plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
	plt.plot(Qprofiles[nCase,:]/meanU[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
	plt.ylim([0,2])
plt.savefig(plume.figdir + 'normProfsI.pdf')
# plt.show()
plt.close()


plt.figure(figsize=(12,6))
clr = plt.cm.PiYG_r(plt.Normalize(-200,200)(plume.read_tag('S',plume.tag)))
clr[..., -1] = plt.Normalize()(fireLine[:,0])
for nCase,Case in enumerate(plume.tag):
	plt.subplot(1,2,2)
	plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
	plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile ')
	# plt.colorbar(label='fireline intensity ')
	plt.ylim([0,2])
	plt.subplot(1,2,1)
	plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
	plt.plot(Qprofiles[nCase,:]/meanU[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
	plt.ylim([0,2])
plt.savefig(plume.figdir + 'normProfsHFX.pdf')
# plt.show()
plt.close()
