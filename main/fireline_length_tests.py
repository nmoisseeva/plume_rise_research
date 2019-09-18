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
testCases = ['W5F7R2','W5F7R2L410km']
#=================end of input===============

runCnt = len(testCases)


plume_tops = np.empty(runCnt) * np.nan
cumT = np.empty(runCnt)* np.nan
plumeTilt = np.empty(runCnt)* np.nan
fireWidth = np.empty((runCnt))* np.nan
parcelHeat = np.empty(runCnt)* np.nan
T0 = np.empty((runCnt,len(plume.lvl))) * np.nan
PMprofiles = np.empty((runCnt,len(plume.lvl)))* np.nan

BLdict = {'Ua':np.empty(runCnt) * np.nan, \
            'Ti':np.empty(runCnt) * np.nan,\
            'zi': np.empty(runCnt)* np.nan}

for nCase,Case in enumerate(testCases):
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
    sliceLoc = wmax_idx[-1]
    PMslicewtop = avedict['pm25'][:,sliceLoc]
    PMprofiles[nCase,:] = PMslicewtop

    #look at what happens with temperature structure----------------------
    #load initial temperature profileT
    if preIgnT:
        T0[nCase] = np.load(profpath)
    else:
        T0[nCase] = avedict['temp'][:,0]

    #BL characteristics -------------------------------------------------
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

    #plot contours
    fig = plt.figure(figsize=(12,6))
    plt.suptitle('%s' %Case)
    plt.subplot(1,2,1)
    plt.title('Time-averaged W')
    # ---w contours and colorbar
    ax1=plt.gca()
    im = ax1.imshow(avedict['w'], origin='lower',extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
    ax1.set_aspect('auto')
    cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    cbari.set_label('horizontal velocity w $[m s^{-2}]$')
    ax1.set_xlabel('horizontal distance [m]')
    ax1.set_ylabel('height AGL [m]')
    ax1.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    ax1.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax1.axhline(y=plume_tops[nCase],ls='--', c='black',label='derived plume top')
    ax1.axvline(x = sliceLoc*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax1.legend()
    # ---non-filled pm contours and colorbar
    cntr = ax1.contour(avedict['pm25'], extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.Greys,levels=pmLevels,linewidths=1)
    # ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
    # cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
    # cbar.set_label('PM2.5 mixing ratio $[ug/kg]$',size=8)
    # cbar.ax.tick_params(labelsize=8)
    # ---heat flux
    axh1 = ax1.twinx()
    axh1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh1.set_ylim([0,150])
    axh1.tick_params(axis='y', colors='red')
    ln = axh1.plot(np.arange(dimX) * plume.dx, avedict['ghfx'], 'r-')
    axh1.set_xlim([0,dimX*plume.dx])


    plt.subplot(1,2,2)
    plt.title('Time-averaged U')
    ax2 = plt.gca()
    # ---u contours and colorbar
    im = ax2.imshow(avedict['u'], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.RdBu_r, vmin=BLdict['Ua'][nCase]-5, vmax=BLdict['Ua'][nCase]+5)
    ax2.set_aspect('auto')
    cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    cbari.set_label('horizontal velocity u $[m s^{-2}]$')
    ax2.set_xlabel('horizontal distance [m]')
    ax2.set_ylabel('height AGL [m]')
    ax2.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    ax2.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax2.axhline(y=plume_tops[nCase],ls='--', c='black',label='derived plume top')
    ax2.axvline(x = sliceLoc*plume.dx, ls=':',c='black',label='location of concentration profile')
    # ---non-filled vapor contours and colorbar
    cntr = ax2.contour(avedict['pm25'], extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.Greys,levels=pmLevels,linewidths=1)
    # ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
    # cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
    # cbar.set_label('PM2.5 mixing ratio $[ug/kg]$',size=8)
    # cbar.ax.tick_params(labelsize=8)
    # ---heat flux
    axh2 = ax2.twinx()
    ln = axh2.plot(np.arange(dimX)*plume.dx, avedict['ghfx'], 'r-')
    axh2.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh2.set_ylim([0,150])
    axh2.tick_params(axis='y', colors='red')
    axh2.set_xlim([0,dimX*plume.dx])

    plt.tight_layout()
    plt.savefig(plume.figdir + 'fireline/ave%s' %Case)
    print('.....-->saved in: %s' %(plume.figdir + 'fireline/ave%s' %Case))
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
print('Regression fit equation: %s ' %regR)
plt.gca().scatter(parcelHeat,cumT)
plt.plot(parcelHeat, regR(parcelHeat))
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.tight_layout()
plt.savefig(plume.figdir + 'fireline/parcelHeat_cumT.pdf')
plt.show()
plt.close()
