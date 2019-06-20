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

preIgnT = 1 		#boolean: whether to use pre-ignition temp profile or averaged upwind profile

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
		avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
	else:
		sys.exit('ERROR: no averaged data found - run prep_plumes.py first!')

	#extract lcoations of max q, w, temp ------------------------------
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
	t_at_wmax = np.array([avedict['temp'][ni,i] for ni, i in enumerate(wmax_idx)])

	#get plume tilt through linear regression of max concentration values
	tilt = np.poly1d(np.polyfit(qmax_idx,plume.lvl[np.isfinite(qmax_profile)],1))

	#create a vertical slize based at a single downwind location of maximum plume rise
	qslicewtop = avedict['qvapor'][:,wmax_idx[-1]] #This should be the proper definition
	Qprofiles[nCase,:] = qslicewtop

	#look at what happens with temperature structure----------------------
	#load initial temperature profileT
	if preIgnT:
		profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
		T0 = np.load(profpath)
	else:
		T0 = avedict['temp'][:,0]

	# #calculate temeperature excess(option1)
	# excessT = t_at_wmax-T0[:len(wmax_idx)]
	# exT[nCase] = np.sum(excessT[excessT>0]) * plume.dz

	# calculate temeprature excess (option 2)
	excessT = tmax_profile[:len(wmax_idx)]-T0[:len(wmax_idx)]
	exT[nCase] = np.sum(excessT[excessT>0]) * plume.dz

	#fire behavior ------------------------------------------------------

	#get fireline intensity and characteristic wind
	# U = np.nanmean(avedict['u'][3:,0]) 	#exclude bottom three layers due to friction  #based on mean wind
	U = avedict['u'][2,0] 					#based on near surface wind
	#store data
	charU[nCase] = U
	print('Wind speed near ground: %s ' %U)

	# #calculate total flux from ground and fire
	ignited = np.array([i for i in avedict['ghfx'] if i > 2])

	#(option 1): just the fire (converted to kinematic)
	totI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 ) #convert to kinematic heat flux (and remove -kilo)

	#(option 2): fire and ground heat flux combined
	# fireI = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005 )
	# grndI = plume.read_tag('S',[Case])[0]* plume.dx * len(ignited)/(1.2 * 1005)
	# totI = fireI + grndI

	#store data
	plumeTilt[nCase] = tilt.c[0]
	fireLine[nCase,:] = totI,len(ignited)
	plume_tops[nCase] = plume.lvl[len(wmax_idx)-1]
	print('plume top: %s' %plume.lvl[len(wmax_idx)-1])



	#estimate atmospheric heating and "cumT --------------------------------"
	# #get "cumulutive" temperature for the profile
	# cumT[nCase] = np.sum(T0[:len(wmax_idx)]*plume.dz)

	#get cumulative T based on  delT
	delT = T0[1:len(wmax_idx)]-T0[0:len(wmax_idx)-1]
	cumT[nCase] = np.sum(delT *	 plume.dz)

	print('Excess temp area along wmax: %d.02; NormI: %d.02' %(exT[nCase],fireLine[nCase,0]/(charU[nCase])))


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

#---------------------calculations over all plumes--------------
normI = fireLine[:,0]/(charU) #normalized fire intensity
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)

print('Variability of normalized intensity for given U level:')
print('%d.2' %np.std(normI))

print('Variability of temperature excess:')
varExcess = np.mean(exT-normI)
print('%d.2 deg' %varExcess)

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
regR0 = np.poly1d(np.polyfit(normI[Rtag==0],cumT[Rtag==0],1))
regR1 = np.poly1d(np.polyfit(normI[Rtag==1],cumT[Rtag==1],1))
# mrkr = ['s' if i==1 else 'o' for i in plume.read_tag('R',RunList)]
print(regR0)
print(regR1)

plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
sc0 = ax.scatter(normI[Rtag==0],cumT[Rtag==0],marker='o', c=plume.read_tag('S',RunList)[Rtag==0], cmap=plt.cm.PiYG_r, vmin=-600, vmax=600)
sc1 = ax.scatter(normI[Rtag==1],cumT[Rtag==1],marker='s', c=plume.read_tag('S',RunList)[Rtag==1], cmap=plt.cm.PiYG_r, vmin=-600, vmax=600)

plt.plot(normI[Rtag==0], regR0(normI[Rtag==0]))
plt.plot(normI[Rtag==1], regR1(normI[Rtag==1]))

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
# clr = plt.cm.PiYG_r(plt.Normalize(-200,200)(plume.read_tag('S',RunList))) 	#color by surface flux
clr = plt.cm.viridis(plt.Normalize(2,12)(charU)) 								#color by windspeed

clr[..., -1] = plt.Normalize()(fireLine[:,0]) 									#add opacity based on fire intensity
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
