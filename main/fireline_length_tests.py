# nmoisseeva@eoas.ubc.ca
# September 2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
from scipy import interpolate
import os.path
import sys
import imp
from matplotlib import ticker


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

preIgnT = 1 		#boolean: whether to use pre-ignition temp profile or averaged upwind profile
plume.fireline_runs = ['W4F7R4L1','W4F7R4','W4F7R4L4']
BLtestruns = ['W4F7R0','W4F7R1','W4F7R3']
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
plt.savefig(plume.figdir + 'fireline/DownwindSmokeProfiles.jpg')
plt.show()

#----------------------------BL TESTS-----------------------------
Up = np.empty((len(BLtestruns),dimZ,dimX)) * np.nan
smokeBL = np.empty((len(BLtestruns),dimZ,dimX)) * np.nan
ziBL = np.empty(len(BLtestruns))
Wp = np.empty((len(BLtestruns),dimZ,dimX)) * np.nan

for nCase,Case in enumerate(BLtestruns):

    #----------check for interpolated data----------------------------
    avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'

    if os.path.isfile(avepath):
        print('Averaged data found at: %s' %avepath)
        avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
    else:
        sys.exit('ERROR: no averaged data found - run prep_plumes.py via submit_interp.sh first on Cedar!')

    U0BL = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')
    T0BL = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')

    #save wind anomaly crossection --------------------------------------
    Uprime = (avedict['u'].T-U0BL).T
    Up[nCase,:,:] = Uprime

    Wp[nCase,:,:] = avedict['w']

    smokeBL[nCase,:,:] = avedict['pm25']            #save smoke

    #save BL height
    dT = T0BL[1:]-T0BL[0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 2
    ziBL[nCase] = plume.dz * (np.argmax(gradT[si:]) + si)


#plot difference between R0 and R1 (same zi, different Gamma)

fig = plt.figure(figsize=(12, 10))
plt.suptitle('INVERSION EFFECT: STRONG (R1) - WEAK (R0)')

plt.subplot(311)
ax1 = plt.gca()
plt.title('U')
im = ax1.imshow(Up[1,:,:]-Up[0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.PiYG_r,vmin=-2, vmax=2)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta u$ $[m/s]$')
ax1.axhline(y = ziBL[0], ls=':', c='darkgrey', label='BL height at ignition')
ax1.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')
ax1.legend()

plt.subplot(312)
ax2 = plt.gca()
plt.title('W')
im = ax2.imshow(Wp[1,:,:]-Wp[0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.PiYG_r,vmin=-2, vmax=2)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta w$ $[m/s]$')
ax2.axhline(y = ziBL[0], ls=':', c='darkgrey', label='BL height at ignition')
ax2.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')

plt.subplot(313)
ax3 = plt.gca()
plt.title('SMOKE')
im = ax3.imshow(smokeBL[1,:,:]-smokeBL[0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.BrBG_r,vmin=-40000, vmax=40000)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta smoke$')
ax3.axhline(y = ziBL[0], ls=':', c='darkgrey', label='BL height at ignition')
ax3.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'fireline/diffGamma.pdf')
plt.savefig(plume.figdir + 'fireline/diffGamma.jpg')

plt.close()



#normalized effects of the zi

fig = plt.figure(figsize=(12, 10))
plt.suptitle('Zi EFFECT: DEEP (R3) - SHALLOW (R1)')

plt.subplot(311)
ax1 = plt.gca()
plt.title('U')
im = ax1.imshow(Up[2,:,:]-Up[1,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.PiYG_r,vmin=-2, vmax=2)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta u$ $[m/s]$')
ax1.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')
ax1.axhline(y = ziBL[2], ls=':', c='black', label='BL height R3')
ax1.axhline(y = ziBL[1], ls='--', c='darkgrey', label='BL height R1')
ax1.legend()

plt.subplot(312)
ax2 = plt.gca()
plt.title('W')
im = ax2.imshow(Wp[2,:,:]-Wp[1,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.PiYG_r,vmin=-2, vmax=2)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta w$ $[m/s]$')
ax2.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')
ax2.axhline(y = ziBL[2], ls=':', c='black', label='BL height R3')
ax2.axhline(y = ziBL[1], ls='--', c='darkgrey', label='BL height R1')

plt.subplot(313)
ax3 = plt.gca()
plt.title('SMOKE')
im = ax3.imshow(smokeBL[2,:,:]-smokeBL[1,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.BrBG_r,vmin=-40000, vmax=40000)
cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
cbari.set_label('$\Delta smoke$')
ax3.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',ylabel='height AGL [m]')
ax3.axhline(y = ziBL[2], ls=':', c='black', label='BL height R3')
ax3.axhline(y = ziBL[1], ls='--', c='darkgrey', label='BL height R1')
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.show()
plt.savefig(plume.figdir + 'fireline/diffZi.pdf')
plt.savefig(plume.figdir + 'fireline/diffZi.jpg')

plt.close()

#----------------------------ROSIE COMPARISON-----------------------------
from scipy.io import netcdf

avepathWRF = plume.wrfdir + 'interp/wrfave_W4F7R4L1.npy'
wrfavedict =  np.load(avepathWRF,allow_pickle=True).item()




pathDALES = plume.wrfdir + 'rosie/Rd00441_20min.nc'
dalesdata = netcdf.netcdf_file(pathDALES, mode ='r')

tWinWRF = 20*15 #seconds to skip for averaging
numWinDALES = int(tWinWRF / 10)

smokeDALES = np.copy(dalesdata.variables['sv_yint'][:,:,:])
aveSVdales = np.mean(smokeDALES[numWinDALES:,:,:],0)

ignited = np.array([i for i in wrfavedict['ghfx'] if i > 2])

#------------Profile Comparison-----------
plt.figure(figsize=(12,4))
plt.subplot(131)
plt.suptitle('TIME-AVERAGED CONCENTRATIONS')
plt.title('1 KM DOWNWIND')
plt.plot(wrfavedict['pm25'][:,50]/np.max(wrfavedict['pm25'][:,50]),plume.lvl)
plt.plot(aveSVdales[:,250]/np.max(aveSVdales[:,250]),dalesdata.variables['zt'][:])
plt.gca().set(xlabel='normalized concentration', ylabel='height AGL [m]')

plt.subplot(132)
plt.title('2 KM DOWNWIND')
plt.plot(wrfavedict['pm25'][:,75]/np.max(wrfavedict['pm25'][:,75]),plume.lvl)
plt.plot(aveSVdales[:,350]/np.max(aveSVdales[:,350]),dalesdata.variables['zt'][:])
plt.gca().set(xlabel='normalized concentration', ylabel='height AGL [m]')
plt.subplot(133)
plt.title('3 KM DOWNWIND')
plt.plot(wrfavedict['pm25'][:,100]/np.max(wrfavedict['pm25'][:,100]),plume.lvl, label='WRF')
plt.plot(aveSVdales[:,450]/np.max(aveSVdales[:,450]),dalesdata.variables['zt'][:], label='DALES')
plt.gca().set(xlabel='normalized concentration', ylabel='height AGL [m]')
plt.legend()
plt.tight_layout(rect=[0, 0, 1, 0.95])

plt.savefig(plume.figdir + 'fireline/ComparisonDALES-WRF_smoke.pdf')
plt.savefig(plume.figdir + 'fireline/ComparisonDALES-WRF_smoke.jpg')

plt.show()
