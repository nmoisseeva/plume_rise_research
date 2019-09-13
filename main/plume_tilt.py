# nmoisseeva@eoas.ubc.ca
# September 2019

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

print('TRACER-BASED ANALYSIS OF LES PLUME DATA')
print('===================================')


param_dict = {'fire':[]}
profile_dict = {'wmax':[],'pm_wmax':[],'tilt':[]}
profile_dict['meta'] = 'wmax: profiles of maximum velocity; \
    pm_wmax: tracer profile at downwind location of maximum vertical velocity \
    tilt: grid number of wmax profile locations'

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
runCnt = len(RunList)

plume_tops = np.empty(runCnt) * np.nan
cumT = np.empty(runCnt)* np.nan
inv_cumT = np.empty(runCnt) * np.nan
plumeTilt = np.empty(runCnt)* np.nan
fireWidth = np.empty((runCnt))* np.nan
parcelHeat = np.empty(runCnt)* np.nan
T0 = np.empty((runCnt,len(plume.lvl))) * np.nan
PMprofiles = np.empty((runCnt,len(plume.lvl)))* np.nan
# exT = np.empty(runCnt)* np.nan
A = np.empty(runCnt)* np.nan

BLdict = {'Ua':np.empty(runCnt) * np.nan, \
            'Ti':np.empty(runCnt) * np.nan,\
            'zi': np.empty(runCnt)* np.nan}

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
        sys.exit('ERROR: no averaged data found - run prep_plumes.py via submit_interp.sh first on Cedar!')

    #extract lcoations of max pm, w, temp ------------------------------
    PMmax_profile = np.nanmax(avedict['pm25'],1) 	#get max q profile
    top_threshold = max(PMmax_profile)*0.001 	#variable threshold (based on near-surface high concentrations!!!!)
    PMmax_profile[PMmax_profile<top_threshold] = np.nan
    PMmax_idx = np.nanargmax(avedict['pm25'][np.isfinite(PMmax_profile)],1)		#get donwind location
    PMmax_meters = PMmax_idx*plume.dx

    wave_plume = avedict['w'].copy()
    wave_plume[avedict['pm25']<top_threshold] = np.nan 		#mask where there is no plume
    wmax_profile = np.nanmax(wave_plume,1) 		#get the profiles
    wmax_idx = np.nanargmax(wave_plume[np.isfinite(wmax_profile)],1)		#get downwind location (index)
    watPM_profile = np.array([avedict['w'][ni,i] for ni, i in enumerate(PMmax_idx)])	#get the profiles

    tmax_profile = np.nanmax(avedict['temp'],1)
    tmax_profile[np.isnan(PMmax_profile)] = np.nan
    tmax_idx = np.nanargmax(avedict['temp'][np.isfinite(tmax_profile)],1)
    watt_profile = np.array([avedict['w'][ni,i] for ni, i in enumerate(tmax_idx)])
    t_at_wmax = np.array([avedict['temp'][ni,i] for ni, i in enumerate(wmax_idx)])

    #get plume tilt through linear regression of max concentration values
    tilt = np.poly1d(np.polyfit(PMmax_idx,plume.lvl[np.isfinite(PMmax_profile)],1))

    #create a vertical slize based at a single downwind location of maximum plume rise
    PMslicewtop = avedict['pm25'][:,wmax_idx[-1]]
    PMprofiles[nCase,:] = PMslicewtop

    #look at what happens with temperature structure----------------------
    #load initial temperature profileT
    if preIgnT:
        profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
        T0[nCase] = np.load(profpath)
    else:
        T0[nCase] = avedict['temp'][:,0]

    # #calculate temeperature excess(option1)
    # excessT = t_at_wmax-T0[:len(wmax_idx)]
    # exT[nCase] = np.sum(excessT[excessT>0]) * plume.dz

    # # calculate temeprature excess (option 2)
    # excessT = tmax_profile[:len(wmax_idx)]-T0[nCase][:len(wmax_idx)]
    # exT[nCase] = np.sum(excessT[excessT>0]) * plume.dz

    #BL characteristics -------------------------------------------------
    # U = np.nanmean(avedict['u'][3:,0]) 	#exclude bottom three layers due to friction  #based on mean wind
    BLdict['Ua'][nCase] = avedict['u'][2,0]
    print('Wind speed near ground: %.1f m/s' %(BLdict['Ua'][nCase]))

    dT = T0[nCase][1:]-T0[nCase][0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 2
    BLdict['Ti'][nCase] = T0[nCase][si+1]              #characteristic BL temperature
    BLdict['zi'][nCase] = plume.dz * (np.argmax(gradT[si:]) + si)

    #fire behavior ------------------------------------------------------

    # #calculate total flux from ground and fire
    ignited = np.array([i for i in avedict['ghfx'] if i > 2])

    #convert to kinematic heat flux (from kW/m2)
    kinI = ignited * 1000/(1.2 * 1005)
    parcelHeat[nCase] = sum(kinI * plume.dx) / BLdict['Ua'][nCase]

    #store data
    plumeTilt[nCase] = tilt.c[0]
    fireWidth[nCase] = len(ignited)
    plume_tops[nCase] = plume.lvl[len(wmax_idx)-1]
    print('Plume top: %d m' %plume_tops[nCase])


    #estimate atmospheric heating and cumT --------------------------------
    #get cumulative T based on  delT
    delTplume = T0[nCase][int(si+1):len(wmax_idx)]-T0[nCase][si:len(wmax_idx)-1]
    cumT[nCase] = np.sum(delTplume * plume.dz)

    # #get invCumT
    # delTfull = T0[nCase][5:]-T0[nCase][4:-1]
    # inv_cumT[nCase] = np.sum(delTfull * plume.dz)

    #charageteristic vertical velocity--------------------------------------------
    g = 9.81
    Fflux = np.mean(ignited) * 1000 / ( 1.2 * 1005)         #get heat flux
    Wf = (g*BLdict['zi'][nCase]*Fflux/BLdict['Ti'][nCase])**(1./3)                     #characteristic velocity based on Deodorff

    # A[nCase] = g* parcelHeat[nCase]/ (zi*Ti[nCase])           #some sort of non-dimensional variable
    #===========================plotting===========================
    #vertical concentration slice at donwind locations of wmax and qmax
    plt.figure(figsize=(12,12))
    plt.suptitle('%s' %Case)
    plt.subplot(2,2,1)
    plt.title('PROFILES OF VERTICAL VELOCITY ALONG W and PM MAXIMA')
    plt.plot(wmax_profile,plume.lvl,'.-',label='$w_{max}$')
    plt.plot(watPM_profile,plume.lvl[np.isfinite(PMmax_profile)],'k.-',label='$w_{qmax}$')
    # plt.axvline(x = Wf, ls='-.', c ='red', label='Wf (fire characteristic velocity)')
    plt.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    plt.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    plt.axhline(y=plume_tops[nCase],ls='--', c='red',label='derived plume top')
    # plt.plot(watt_profile,plume.lvl[np.isfinite(tmax_profile)],'r.-',label='$w_{tmax}$')
    plt.xlabel('velocity [m/s]')
    plt.ylabel('height [m]')
    plt.legend()
    plt.ylim([0,plume.lvl[-1]])

    plt.subplot(2,2,2)
    plt.title('HORIZONTAL LOCATION OF EXTREMA')
    plt.plot(wmax_idx,plume.lvl[np.isfinite(wmax_profile)],'.-',label='$w_{max}$')
    plt.plot(PMmax_idx,plume.lvl[np.isfinite(PMmax_profile)],'k.--',label='$q_{max}$')
    # plt.plot(tmax_idx,plume.lvl[np.isfinite(tmax_profile)],'r.--',label='$t_{max}$')
    plt.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    plt.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    plt.axhline(y=plume_tops[nCase],ls='--', c='red',label='derived plume top')
    plt.plot(PMmax_idx, tilt(PMmax_idx))
    plt.xlabel('x distance [m]')
    ax = plt.gca()
    ax.set_xticks(np.arange(0,120,20))
    ax.set_xticklabels(np.arange(0,120,20)*40)
    plt.ylabel('height [m]')
    plt.ylim([0,plume.lvl[-1]])
    plt.legend()

    plt.subplot(2,2,3)
    plt.title('PM CONCENTRATION DOWNWIND (SLICE)')
    plt.plot(PMslicewtop, plume.lvl, 'g.--', label='concentration based on $w_{top}$')
    plt.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    plt.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    plt.axhline(y=plume_tops[nCase],ls='--', c='red',label='derived plume top')
    plt.xlabel('PM2.5 concentration [ug/kg]')
    plt.ylabel('height [m]')
    plt.ylim([0,plume.lvl[-1]])
    plt.legend()

    plt.subplot(2,2,4)
    plt.title('PLUME vs AMBIENT TEMPERATURE')
    plt.plot(T0[nCase], plume.lvl, label='pre-ignition profile',c='lightblue')
    plt.plot(t_at_wmax,plume.lvl[:len(wmax_idx)],c = 'orange',label='in-plume temperature')
    plt.plot(avedict['temp'][:len(wmax_idx),0],plume.lvl[:len(wmax_idx)],label='ambient temperature')
    plt.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    plt.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    plt.axhline(y=plume_tops[nCase],ls='--', c='red',label='derived plume top')
    plt.xlabel('temperature [K]')
    plt.ylabel('height [m]')
    plt.legend()
    plt.ylim([0,plume.lvl[-1]])
    plt.xlim([285,330])
    # plt.show()
    plt.savefig(plume.figdir + 'profiles_%s.pdf' %Case)
    plt.close()

#---------------------calculations over all plumes--------------
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)

#WHAT LEVEL TO PICK FOR defining characteristic U
# print('Variability of normalized intensity for given U level:')
# print('%d.2' %np.std(parcelHeat))

#-------------subplots of tilt predictors-----------
#create scatter plots of tilts vs other variables
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('FIRELINE INTENSITY vs TILT')
ax1 = plt.scatter(parcelHeat,plumeTilt,c=BLdict['Ua'],cmap=plt.cm.viridis)
plt.xlabel('fireline intensity [kW/m]')
plt.ylabel('plume tilt')
plt.colorbar(ax1,label='windspeed [m/s]')

plt.subplot(1,3,2)
plt.title('WIND vs TILT')
ax2 = plt.scatter(BLdict['Ua'],plumeTilt,c=parcelHeat,cmap=plt.cm.plasma)
plt.xlabel('mean wind (m/s)')
plt.ylabel('plume tilt')
plt.colorbar(ax2, label='fireline intensity [kW/m]')

plt.subplot(1,3,3)
plt.title('WIDTH vs TILT')
ax3 = plt.scatter(fireWidth,plumeTilt,c=BLdict['Ua'],cmap=plt.cm.viridis)
plt.xlabel('fireline width [#grids]')
plt.ylabel('plume tilt')
plt.colorbar(ax3,label='windspeed [m/s]')

plt.tight_layout()
plt.savefig(plume.figdir + 'tilt_subplots.pdf')
# plt.show()
plt.close()


#-------------subplots for eatch R run-----------
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('R0 RUNS')
ax1 = plt.scatter(parcelHeat[Rtag==0],cumT[Rtag==0],marker='o', c=plume.read_tag('W',RunList)[Rtag==0], cmap=plt.cm.viridis)
for i, txt in enumerate(plume.read_tag('F',RunList)[Rtag==0]):
    plt.gca().annotate(txt, (parcelHeat[Rtag==0][i]+20,cumT[Rtag==0][i]+20), fontsize=9)
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.colorbar(ax1,label='windspeed [m/s]')
plt.xlim([0,6000])
plt.ylim([0,1000])

plt.subplot(1,3,2)
plt.title('R1 RUNS')
ax2 = plt.scatter(parcelHeat[Rtag==1],cumT[Rtag==1],marker='s', c=plume.read_tag('W',RunList)[Rtag==1], cmap=plt.cm.viridis)
for i, txt in enumerate(plume.read_tag('F',RunList)[Rtag==1]):
    plt.gca().annotate(txt, (parcelHeat[Rtag==1][i]+20,cumT[Rtag==1][i]+20), fontsize=9)
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.colorbar(ax1,label='windspeed [m/s]')
plt.xlim([0,6000])
plt.ylim([0,1000])

plt.subplot(1,3,3)
plt.title('R2 RUNS')
ax3 = plt.scatter(parcelHeat[Rtag==2],cumT[Rtag==2],marker='^', c=plume.read_tag('W',RunList)[Rtag==2], cmap=plt.cm.viridis)
for i, txt in enumerate(plume.read_tag('F',RunList)[Rtag==2]):
    plt.gca().annotate(txt, (parcelHeat[Rtag==2][i]+20,cumT[Rtag==2][i]+20), fontsize=9)
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.colorbar(ax1,label='windspeed [m/s]')
plt.xlim([0,6000])
plt.ylim([0,1000])

plt.tight_layout()
plt.savefig(plume.figdir + 'Rn_sublots.pdf')
plt.show()
plt.close()

#-------------main plot with CumT-----------
#get linear regression for the normalized fireLine
regR0 = np.poly1d(np.polyfit(parcelHeat[Rtag==0],cumT[Rtag==0],1))
regR1 = np.poly1d(np.polyfit(parcelHeat[Rtag==1],cumT[Rtag==1],1))
regR2 = np.poly1d(np.polyfit(parcelHeat[Rtag==2],cumT[Rtag==2],1))
# mrkr = ['s' if i==1 else 'o' for i in plume.read_tag('R',RunList)]
print('R0 profiles: regression fit equation: %s ' %regR0)
print('R1 profiles: regression fit equation: %s ' %regR1)
print('R2 profiles: regression fit equation: %s ' %regR2)

plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
sc0 = ax.scatter(parcelHeat[Rtag==0],cumT[Rtag==0],marker='o', c='C0',label='R0')
sc1 = ax.scatter(parcelHeat[Rtag==1],cumT[Rtag==1],marker='s', c='C1',label='R1')
sc2 = ax.scatter(parcelHeat[Rtag==2],cumT[Rtag==2],marker='^', c='C2',label='R2')

plt.plot(parcelHeat[Rtag==0], regR0(parcelHeat[Rtag==0]), 'C0')
plt.plot(parcelHeat[Rtag==1], regR1(parcelHeat[Rtag==1]), 'C1')
plt.plot(parcelHeat[Rtag==2], regR2(parcelHeat[Rtag==2]), 'C2')
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.legend()
plt.tight_layout()
plt.savefig(plume.figdir + 'parcelHeat_cumT.pdf')
plt.show()
plt.close()



#
#
# plt.figure(figsize=(12,6))
# clr = plt.cm.plasma(plt.Normalize()(parcelHeat))
# # clr[..., -1] = plt.Normalize()(fireLine[:,0])
# for nCase,Case in enumerate(RunList):
#     plt.subplot(1,2,2)
#     plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
#     plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile ')
#     # plt.colorbar(label='fireline intensity ')
#     plt.ylim([0,2])
#     plt.subplot(1,2,1)
#     plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
#     plt.plot(Qprofiles[nCase,:]/Ua[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
#     plt.ylim([0,2])
# plt.savefig(plume.figdir + 'normProfsI.pdf')
# # plt.show()
# plt.close()
#
#
# plt.figure(figsize=(12,6))
# # clr = plt.cm.PiYG_r(plt.Normalize(-200,200)(plume.read_tag('S',RunList))) 	#color by surface flux
# clr = plt.cm.viridis(plt.Normalize(2,12)(Ua)) 								#color by windspeed
#
# clr[..., -1] = plt.Normalize()(parcelHeat) 									#add opacity based on fire intensity
# for nCase,Case in enumerate(RunList):
#     plt.subplot(1,2,2)
#     plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
#     plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_tops[nCase],c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile ')
#     # plt.colorbar(label='fireline intensity ')
#     plt.ylim([0,2])
#     plt.subplot(1,2,1)
#     plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
#     plt.plot(Qprofiles[nCase,:]/Ua[nCase],plume.lvl/plume_tops[nCase], c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
#     plt.ylim([0,2])
# plt.savefig(plume.figdir + 'normProfsHFX.pdf')
# # plt.show()
# plt.close()

plt.title('INITAL ATMOSPHERIC PROFILES (R0 | R1 | R2)')
for Case in T0[Rtag==0]:
    lR0 = plt.plot(Case, plume.lvl, color='C0', label='R0')
for Case in T0[Rtag==1]:
    lR1 = plt.plot(Case, plume.lvl, color ='C1', label='R1')
for Case in T0[Rtag==2]:
    lR2 = plt.plot(Case, plume.lvl, color ='C2', label='R2')
plt.xlabel('potential temperature [K]')
plt.ylabel('height [m]')
plt.legend(handles=[lR0[-1],lR1[-1],lR2[-1]])
plt.savefig(plume.figdir + 'T0profiles.pdf')
plt.show()

# plt.scatter(Ti,A)
# plt.show()
