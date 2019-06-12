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

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
runCnt = len(RunList)

plume_tops = np.empty(runCnt) * np.nan
cumT = np.empty(runCnt)* np.nan
charU = np.empty(runCnt) * np.nan		#characteristic wind speed (could be mean BL or near surface)
plumeTilt = np.empty(runCnt)* np.nan
fireLine = np.empty((runCnt,2))* np.nan #intensity, width
Qprofiles = np.empty((runCnt,len(plume.lvl)))* np.nan
exT = np.empty(runCnt)* np.nan

for nCase,Case in enumerate(RunList):
	if Case in plume.exclude_runs:
		continue
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
	top_threshold = max(qmax_profile)*0.001 	#variable threshold
	# top_threshold = 0.1 						#hardcoded threshold
	qmax_profile[qmax_profile<top_threshold] = np.nan
	qmax_idx = np.nanargmax(avedict['qvapor'][np.isfinite(qmax_profile)],1)		#get donwind location
	qmax_meters = qmax_idx*plume.dx

	wave_plume = avedict['w'].copy()
	wave_plume[avedict['qvapor']<top_threshold] = np.nan 		#mask where there is no plume
	wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
	wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
	watq_profile = np.array([avedict['w'][ni,i] for ni, i in enumerate(qmax_idx)])	#get the profiles

	tmax_profile = np.nanmax(avedict['temp'],1)
	tmax_profile[np.isnan(qmax_profile)] = np.nan
	tmax_idx = np.nanargmax(avedict['temp'][np.isfinite(tmax_profile)],1)
	watt_profile = np.array([avedict['w'][ni,i] for ni, i in enumerate(tmax_idx)])

	#calculate temeperature excess(option1)
	t_at_wmax = np.array([avedict['temp'][ni,i] for ni, i in enumerate(wmax_idx)])
	# #check if max temp location near ground is same as max velocity - it's not most of the time!!!! but this desn't mean we should use it
	# surfaceT_idx = np.nanargmax(avedict['temp'][0,:])
	# if surfaceT_idx != wmax_idx[0]:
	# 	t_at_wmax[0] = avedict['temp'][0,surfaceT_idx]
	exT[nCase] = np.sum(t_at_wmax-avedict['temp'][:len(wmax_idx),0])
	# #calculate temeprature excess (option 2)
	# t_at_wmax = np.array([avedict['temp'][ni,i] for ni, i in enumerate(wmax_idx)])
	# tmax_profile = np.nanmax(avedict['temp'],1)
	# exT[nCase] = np.sum(tmax_profile[:len(wmax_idx)]-avedict['temp'][:len(wmax_idx),0])



	# #regF values
	# 0.1 = poly1d([1.27427632e+02, 3.30092886e+05])
	#0.01 = 124.4 x + 3.774e+05
	#0.2 = 130.8 x + 3.051e+05

	#get plume tilt through linear regression of max concentration values
	tilt = np.poly1d(np.polyfit(qmax_idx,plume.lvl[np.isfinite(qmax_profile)],1))

	#get fireline intensity and wind
	# U = np.nanmean(avedict['u'][3:,0]) 	#exclude bottom three layers due to friction  #based on mean wind
	U = avedict['u'][2,0] 					#based on near surface wind
	#store data
	charU[nCase] = U
	print('Wind speed near ground: %s ' %U)

	ignited = np.array([i for i in avedict['ghfx'] if i > 2])
	# #calculate total flux from ground and fire
	# fireI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 )
	# grndI = plume.read_tag('S',[Case])[0]* plume.dx * len(ignited)/(1.2 * 1005)
	# totI = fireI + grndI
	totI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 ) #convert to kinematic heat flux (and remove -kilo)
	# totI = sum(ignited) * 1000 / ( 1.2 * 1005 ) #convert to kinematic heat flux (and remove -kilo)


	plumeTilt[nCase] = tilt.c[0]
	fireLine[nCase,:] = totI,len(ignited)
	plume_tops[nCase] = plume.lvl[len(wmax_idx)-1]
	print('plume top: %s' %plume.lvl[len(wmax_idx)-1])


	#create a vertical slize based at a single downwind location of maximum plume rise
	qslicewtop = avedict['qvapor'][:,wmax_idx[-1]] #This should be the proper definition

	#get "cumulutive" temperature for the profile
	cumT[nCase] = np.sum(avedict['temp'][:len(wmax_idx),0]*plume.dz)
	#
	# #get cumulative T based on  delT
	# delT = avedict['temp'][1:len(wmax_idx),0]-avedict['temp'][0:len(wmax_idx)-1,0]
	# cumT[nCase] = np.sum(delT *	 plume.dz)

	# #get cumulative T based on  delT
	# delT = avedict['temp'][1:len(wmax_idx),0]-avedict['temp'][0:len(wmax_idx)-1,0]
	# cumT[nCase] = np.sum(delT) + avedict['temp'][0,0]

	# #get "cumulutive" temperature for the profile assuming uniform grid (in both dx and dz)
	# cumT[nCase] = np.sum(avedict['temp'][:len(wmax_idx),0])

	print('Excess temp along wmax: %d.02; NormI/dx: %d.02' %(exT[nCase],fireLine[nCase,0]/(charU[nCase]*plume.dx)))


	#===========================plotting===========================
	#vertical concentration slice at donwind locations of wmax and qmax
	plt.figure(figsize=(12,12))
	plt.suptitle('%s' %Case)
	plt.subplot(2,2,1)
	plt.title('VALUE ALONG MAX PROFILE OF W, MAX AND MIN U')
	plt.plot(wmax_profile,plume.lvl,'.-',label='$w_{max}$')
	plt.plot(watq_profile,plume.lvl[np.isfinite(qmax_profile)],'k.-',label='$w_{qmax}$')
	# plt.plot(watt_profile,plume.lvl[np.isfinite(tmax_profile)],'r.-',label='$w_{tmax}$')
	plt.xlabel('velocity [m/s]')
	plt.ylabel('height [m]')
	plt.legend()
	plt.ylim([0,plume.lvl[-1]])

	plt.subplot(2,2,2)
	plt.title('HORIZONTAL LOCATION OF EXTREMA')
	plt.plot(wmax_idx,plume.lvl[np.isfinite(wmax_profile)],'.-',label='$w_{max}$')
	plt.plot(qmax_idx,plume.lvl[np.isfinite(qmax_profile)],'k.--',label='$q_{max}$')
	plt.plot(tmax_idx,plume.lvl[np.isfinite(tmax_profile)],'r.--',label='$t_{max}$')
	plt.plot(qmax_idx, tilt(qmax_idx))
	plt.xlabel('x distance [m]')
	ax = plt.gca()
	ax.set_xticks(np.arange(0,75,15))
	ax.set_xticklabels(np.arange(0,75,15)*40)
	plt.ylabel('height [m]')
	plt.ylim([0,plume.lvl[-1]])
	plt.legend()

	plt.subplot(2,2,3)
	plt.title('Q CONCENTRATION DOWNWIND (SLICE)')
	plt.plot(qslicewtop, plume.lvl, 'g.--', label='based on $w_{top}$')
	plt.xlabel('water vapor [g/kg]')
	plt.ylabel('height [m]')
	plt.ylim([0,plume.lvl[-1]])
	plt.legend()

	plt.subplot(2,2,4)
	plt.title('PLUME vs AMBIENT TEMPERATURE')
	plt.plot(t_at_wmax,plume.lvl[:len(wmax_idx)],c = 'orange',label='in-plume temperature')
	plt.plot(avedict['temp'][:len(wmax_idx),0],plume.lvl[:len(wmax_idx)],label='ambient temperature')
	plt.xlabel('temperature [K]')
	plt.ylabel('height [m]')
	plt.legend()
	plt.ylim([0,plume.lvl[-1]])
	plt.xlim([285,330])
	# plt.show()
	plt.savefig(plume.figdir + 'profiles_%s.pdf' %Case)
	plt.close()


	Qprofiles[nCase,:] = qslicewtop
	# plt.plot(avedict['temp'][0:len(wmax_idx),0],plume.lvl[:len(wmax_idx)])
	# plt.ylim([0,2500])
# plt.show()

#---------------------calculations over all plumes--------------
normI = fireLine[:,0]/(charU) #normalized fire intensity
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)

print('Variability of normalized intensity:')
print(np.std(normI))

#-------------subplots of tilt predictors-----------
#create scatter plots of tilts vs other variables
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('FIRELINE INTENSITY vs TILT')
ax1 = plt.scatter(fireLine[:,0],plumeTilt,c=charU,cmap=plt.cm.viridis)
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('plume tilt')
plt.colorbar(ax1,label='windspeed [m/s]')

plt.subplot(1,3,2)
plt.title('WIND vs TILT')
ax2 = plt.scatter(charU,plumeTilt,c=fireLine[:,0],cmap=plt.cm.plasma)
plt.xlabel('mean wind (m/s)')
plt.ylabel('plume tilt')
plt.colorbar(ax2, label='fireline intensity [kW/m]')

plt.subplot(1,3,3)
plt.title('WIDTH vs TILT')
ax3 = plt.scatter(fireLine[:,1],plumeTilt,c=charU,cmap=plt.cm.viridis)
plt.xlabel('fireline width [#grids]')
plt.ylabel('plume tilt')
plt.colorbar(ax3,label='windspeed [m/s]')

plt.tight_layout()
plt.savefig(plume.figdir + 'tilt_subplots.pdf')
# plt.show()
plt.close()



plt.title('FIRELINE INTENSITY vs CumT')
ax1 = plt.scatter(fireLine[:,0],cumT, c=charU)
# ax = plt.gca()
# ax.scatter(fireLine[:,0],cumT, c=charU)
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
# mrkr = ['s' if i==1 else 'o' for i in plume.read_tag('R',RunList)]
print(regF)

plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
sc0 = ax.scatter(normI[Rtag==0],cumT[Rtag==0],marker='o', c=plume.read_tag('S',RunList)[Rtag==0], cmap=plt.cm.PiYG_r, vmin=-600, vmax=600)
sc1 = ax.scatter(normI[Rtag==1],cumT[Rtag==1],marker='s', c=plume.read_tag('S',RunList)[Rtag==1], cmap=plt.cm.PiYG_r, vmin=-600, vmax=600)

plt.plot(normI, regF(normI))
for i, txt in enumerate(plume.read_tag('F',RunList)):
    ax.annotate(txt, (normI[i]+100,cumT[i]+100), fontsize=9)
plt.colorbar(sc0, label='surface heat flux [$W/m^{2}$]')
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.tight_layout()
plt.savefig(plume.figdir + 'normI_cumT.pdf')
plt.show()
plt.close()










plt.figure(figsize=(12,6))
clr = plt.cm.plasma(plt.Normalize()(fireLine[:,0]))
# clr[..., -1] = plt.Normalize()(fireLine[:,0])
for nCase,Case in enumerate(RunList):
	plt.subplot(1,2,2)
	plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
	plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile ')
	# plt.colorbar(label='fireline intensity ')
	plt.ylim([0,2])
	plt.subplot(1,2,1)
	plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
	plt.plot(Qprofiles[nCase,:]/charU[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
	plt.ylim([0,2])
plt.savefig(plume.figdir + 'normProfsI.pdf')
# plt.show()
plt.close()


plt.figure(figsize=(12,6))
clr = plt.cm.PiYG_r(plt.Normalize(-200,200)(plume.read_tag('S',RunList)))

clr[..., -1] = plt.Normalize()(fireLine[:,0])
for nCase,Case in enumerate(RunList):
	plt.subplot(1,2,2)
	plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
	plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile ')
	# plt.colorbar(label='fireline intensity ')
	plt.ylim([0,2])
	plt.subplot(1,2,1)
	plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
	plt.plot(Qprofiles[nCase,:]/charU[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
	plt.ylabel('normalized height')
	plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
	plt.ylim([0,2])
plt.savefig(plume.figdir + 'normProfsHFX.pdf')
# plt.show()
plt.close()
