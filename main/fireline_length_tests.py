# nmoisseeva@eoas.ubc.ca
# September 2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import sys
import imp

#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

preIgnT = 1 		#boolean: whether to use pre-ignition temp profile or averaged upwind profile
plume.fireline_runs = ['W4F7R4L1','W4F7R4','W4F7R4L4']
#=================end of input===============

runCnt = len(plume.fireline_runs)


plume_tops = np.empty(runCnt) * np.nan
cumT = np.empty(runCnt)* np.nan
plumeTilt = np.empty(runCnt)* np.nan
fireWidth = np.empty((runCnt))* np.nan
parcelHeat = np.empty(runCnt)* np.nan
T0 = np.empty((runCnt,len(plume.lvl))) * np.nan
U0 = np.empty((runCnt,len(plume.lvl))) * np.nan
PMprofiles = np.empty((runCnt,len(plume.lvl)))* np.nan
UprimeCS=[]
SmokeProfiles=np.empty((runCnt,3,len(plume.lvl)))* np.nan

BLdict = {'Ua':np.empty(runCnt) * np.nan, \
            'Ti':np.empty(runCnt) * np.nan,\
            'zi': np.empty(runCnt)* np.nan}

for nCase,Case in enumerate(plume.fireline_runs):
    #----------get preignition temperature profile-----------------
    cspath = plume.wrfdir + 'interp/wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    #save initial temperature prfile
    profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
    profileT = np.mean(csdict['temp'][0,:,:],1)
    np.save(profpath,profileT)

    #----------check for interpolated data----------------------------
    avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'

    if os.path.isfile(avepath):
        print('Averaged data found at: %s' %avepath)
        avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
    else:
        sys.exit('ERROR: no averaged data found - run prep_plumes.py via submit_interp.sh first on Cedar!')

    #extract lcoations of max pm, w, temp ------------------------------
    PMmax_profile = np.nanmax(avedict['pm25'],1) 	#get max smoke profile
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
    sliceLoc = wmax_idx[-1]
    PMslicewtop = avedict['pm25'][:,sliceLoc]
    PMprofiles[nCase,:] = PMslicewtop

    #look at what happens with temperature structure----------------------
    #load initial profiles
    if preIgnT:
        T0[nCase] = np.load(profpath)
    else:
        T0[nCase] = avedict['temp'][:,0]

    U0[nCase] = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    #BL characteristics -------------------------------------------------
    BLdict['Ua'][nCase] = avedict['u'][2,0]
    print('Wind speed near ground: %.1f m/s' %(BLdict['Ua'][nCase]))

    dT = T0[nCase][1:]-T0[nCase][0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 2
    BLdict['Ti'][nCase] = T0[nCase][si+1]              #characteristic BL temperature
    BLdict['zi'][nCase] = plume.dz * (np.argmax(gradT[si:]) + si)

    #save wind anomaly crossection --------------------------------------
    Uprime = (avedict['u'].T-U0[nCase]).T
    ReverseFlow = Uprime[10,:]
    UprimeCS.append(ReverseFlow)
    print(avedict['ghfx']  > 2)

    #save smoke profiles  --------------------------------------
    SmokeProfiles[nCase,0,:] = avedict['pm25'][:,55]/np.max(avedict['pm25'][:,55])
    SmokeProfiles[nCase,1,:] = avedict['pm25'][:,100]/np.max(avedict['pm25'][:,100])
    SmokeProfiles[nCase,2,:] = avedict['pm25'][:,150]/np.max(avedict['pm25'][:,150])


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

    #===========================plotting===========================

    maxPM = int(np.max(avedict['pm25']))
    pmLevels = np.arange(maxPM*0.001,maxPM/2.,maxPM/10.)
    uLevels = np.arange(int(BLdict['Ua'][nCase]) - 5, int(BLdict['Ua'][nCase])+5, 1)
    dimZ, dimX = np.shape(avedict['w'])
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
    plt.savefig(plume.figdir + 'fireline/profiles_%s.pdf' %Case)
    plt.close()

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
plt.savefig(plume.figdir + 'fireline/tilt_subplots.pdf')
# plt.show()
plt.close()


#-------------main plot with CumT-----------
plt.figure()
plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
regR = np.poly1d(np.polyfit(parcelHeat,cumT,1))
plt.gca().scatter(parcelHeat,cumT)
plt.plot(parcelHeat, regR(parcelHeat))
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
for i, txt in enumerate(plume.fireline_runs):
    plt.gca().annotate(txt, (parcelHeat[i]+2,cumT[i]+2), fontsize=9)
plt.tight_layout()
plt.savefig(plume.figdir + 'fireline/parcelHeat_cumT.pdf')
plt.close()

#------------Velocity Enhancement-----------
plt.figure()
haxis = np.arange(dimX)*plume.dx
ax = plt.gca()
plt.title('FIRE-INDUCED WIND DYNAMCIS at 400 m AGL')
plt.plot(haxis, UprimeCS[0], label='1 km')
plt.plot(haxis, UprimeCS[1], label='2 km')
plt.plot(haxis, UprimeCS[2], label='4 km')
# plt.axhline(y=0,xmin=0,xmax=16000,color='grey')
plt.xlabel('horizontal distance [m]')
plt.ylabel('Uprime [m/s]' )
plt.xlim([0,max(haxis)])
plt.tight_layout()
plt.legend()
plt.savefig(plume.figdir + 'fireline/Uprime400m.pdf')
plt.show()
plt.close()

#------------Profile Comparison-----------
plt.figure(figsize=(12,4))
plt.subplot(131)
plt.suptitle('TIME-AVERAGED CONCENTRATIONS')
plt.title('1 KM DOWNWIND')
plt.plot(SmokeProfiles[0,0,:],plume.lvl)
plt.plot(SmokeProfiles[1,0,:],plume.lvl)
plt.plot(SmokeProfiles[2,0,:],plume.lvl)
plt.gca().set(xlabel='normalized concentration', ylabel='height [m]')
plt.subplot(132)
plt.title('3 KM DOWNWIND')
plt.plot(SmokeProfiles[0,1,:],plume.lvl)
plt.plot(SmokeProfiles[1,1,:],plume.lvl)
plt.plot(SmokeProfiles[2,1,:],plume.lvl)
plt.gca().set(xlabel='normalized concentration', ylabel='height [m]')
plt.subplot(133)
plt.title('5 KM DOWNWIND')
plt.plot(SmokeProfiles[0,2,:],plume.lvl,label='1 km fireline')
plt.plot(SmokeProfiles[1,2,:],plume.lvl,label='2 km fireline')
plt.plot(SmokeProfiles[2,2,:],plume.lvl,label='4 km fireline')
plt.gca().set(xlabel='normalized concentration', ylabel='height [m]')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.legend()
plt.savefig(plume.figdir + 'fireline/DownwindSmokeProfiles.pdf')
plt.show()
