# January 2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.stats import linregress



#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time


#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
# RunList = ['W5F6R3','W4F7R1','W4F7R4L1','W4F7R4','W4F7R4L4']


runCnt = len(RunList)
g = 9.81


#take a vertical cross-section (ensure it's stable)
#calculate CDF of cumT (delT*dz)
#calculate temperature gradient -
#map temperature gradient to CDF

#save variables for dimensional analysis
r = np.empty((runCnt)) * np.nan
Ua = np.empty((runCnt)) * np.nan
zi = np.empty((runCnt)) * np.nan
zCL = np.empty((runCnt)) * np.nan
Omega = np.empty((runCnt)) * np.nan
Phi = np.empty((runCnt)) * np.nan
FI = np.empty((runCnt)) * np.nan
width = np.empty((runCnt)) * np.nan
FirelineProfiles, FirelineUprime400m = [], []

for nCase,Case in enumerate(RunList):
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)

    csdict = plume.load_CS_prep_Profiles(Case)
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    #mask plume with 30ppm---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= 30, csdict['pm25'][-1,:,:] )
    zi[nCase] = plume.get_zi(T0)
    si = 3

    #locate centerline
    ctrZidx = pm.argmax(0)
    ctrXidx = pm.argmax(1)
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])
    for nZ in range(dimZ):
        if pmCtr[ctrXidx[nZ]] < pm[nZ,ctrXidx[nZ]]:
            pmCtr[ctrXidx[nZ]] = pm[nZ,ctrXidx[nZ]]


    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
    centerline = ma.masked_where(plume.lvl[ctrZidx] == 0, plume.lvl[ctrZidx])
    tilt = ymax/xmax                #tilt is calculated based on max centerline height
    smoothCenterline = savgol_filter(centerline, 31, 3) # window size 101, polynomial order 3


    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 101, 3) # window size 101, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.05 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]
    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.mean(stablePM,1)
    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)


    #define source (r and H)------------------------
    #raduis using full 2D average - THIS IS THE APPROACH USED FOR DIMENSIONAL ANALYSIS
    masked_flux = ma.masked_less_equal(csdict['ghfx2D'], 0)
    cs_flux = np.nanmean(masked_flux,1)                         #get cross section for each timestep
    fire = []
    xmax = np.argmax(cs_flux,axis=1)                            #get max for each timestep
    for nP, pt in enumerate(xmax[plume.ign_over:]):             #excludes steps containing ignition
        subset = cs_flux[plume.ign_over+nP,pt-plume.wi:pt+plume.wf]     #set averaging window around a maximum
        fire.append(subset)
    burning = np.sum(np.sum(csdict['ghfx2D'][plume.ign_over:,:,:],1),1) #get total heat flux from entire fire area
    firedepth_t = np.sum(csdict['ghfx2D'][plume.ign_over:,:,:]>0,2)
    maskeddepth = ma.masked_less_equal(firedepth_t,0)
    width[nCase] = np.mean(np.mean(maskeddepth,1))

    meanFire = np.nanmean(fire,0)
    ignited = np.array([i for i in meanFire if i > 0.5])
    H = np.mean(ignited) * 1000 / ( 1.2 * 1005)         #get heat flux
    r[nCase] = len(ignited) * plume.dx
    Phi[nCase] = sum(ignited*plume.dx) * 1000 / ( 1.2 * 1005)
    FI[nCase] = np.mean(burning)*  1000 / ( 1.2 * 1005)

    #compare with center average only
    avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'
    avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
    ignitedCtr = np.array([i for i in avedict['ghfx'] if i > 0.5])
    HCtr = np.mean(ignitedCtr) * 1000 / ( 1.2 * 1005)         #get heat flux
    rCtr = len(ignitedCtr) * plume.dx

    #plot for reference
    whaxis = np.arange(len(meanFire))*plume.dx
    plt.title('SLAB vs. CENTER AVERAGE: %s' %Case)
    ax = plt.gca()
    plt.plot(whaxis[:150],meanFire[:150], label='slab average: r = %s, H = %2d' %(r[nCase],H))
    # ax.fill_between(whaxis[:150], 0, 1, where=meanFire[:150]>0.5, color='red', alpha=0.1, transform=ax.get_xaxis_transform(), label='averaging window')
    plt.plot(whaxis[:150],avedict['ghfx'][:150],color='C1',linestyle='--',label='center average:r = %s, H = %2d' %(rCtr,HCtr))
    # ax.fill_between(whaxis[:150], 0, 1, where=avedict['ghfx'][:150]>0.5, color='grey', alpha=0.1, transform=ax.get_xaxis_transform(), label='averaging window')
    ax.set(ylabel='heat flux [kW/m2]')
    plt.legend()
    plt.savefig(plume.figdir + 'fireDiagnostics/fire%s.pdf' %Case)
    plt.close()
    # plt.show()

    #dimensional analysis variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])

    zCLidx = int(np.mean(ctrZidx[1:][stablePMmask]))
    dT = T0[1:]-T0[0:-1]
    Omega[nCase] = np.sum(dT[si+1:zCLidx]*plume.dz)
    if Omega[nCase] < 0 :
        print('\033[93m' + '$\Omega$: %0.2f ' %Omega[nCase] + '\033[0m')
    Ua[nCase] = np.mean(U0[si:zCLidx])


    #PLOTTING =========================================================
    haxis = np.arange(dimX)*plume.dx
    maxPM = int(np.max(csdict['pm25']))
    pmLevels = np.geomspace(30,maxPM/10,num=10)
    PMcontours = ma.masked_where(csdict['pm25'] <= 30,csdict['pm25'])

    fig = plt.figure(figsize=(22,8))
    gs = fig.add_gridspec(ncols=2, nrows=2,width_ratios=[6,1])
    plt.suptitle('%s' %Case)

    ax1=fig.add_subplot(gs[0])
    axh1=ax1.twinx()
    # ---u contours and colorbar
    im = ax1.imshow((csdict['u'][-1,:,:].T - U0).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)
    cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    cbari.set_label('ralative horizontal velocity $[m s^{-1}]$')
    # ---non-filled vapor contours and colorbar
    cntr = ax1.contour(PMcontours[-1,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
    ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    ax1.axhline(y = zi[nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax1.set(ylabel='height AGL [m]')
    ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
    ax1.legend()
    # ---heat flux
    ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][-1,:], 'r-')
    axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
    axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])
    axh1.tick_params(axis='y', colors='red')


    ax2=fig.add_subplot(gs[1])
    fim = ax2.imshow(csdict['ghfx2D'][-1,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
    cbarif = fig.colorbar(fim, orientation='horizontal')
    cbarif.set_label('heat flux [$kW / m^2$]')
    ax2.set(xlabel='x distance [m]',ylabel='y distance [m]',aspect='equal')


    ax3=fig.add_subplot(gs[2])
    l1, = plt.plot(haxis, pmCtr, label='concentration along centerline', color='C1')
    l2 = ax3.fill_between(haxis[1:], 0, 1, where=stablePMmask, color='grey', alpha=0.4, transform=ax3.get_xaxis_transform(), label='averaging window')
    ax3.set(xlim=[0,dimX*plume.dx],xlabel='horizontal distance [m]',ylabel='concentration [ug/kg]' )
    ax32 = ax3.twinx()
    l3, = plt.plot(haxis,smoothCenterline, label='smoothed centerline height ', color='C2',linestyle=':')
    ax32.set(xlim=[0,dimX*plume.dx],ylim=[0,2000], ylabel='height [m]' )
    plt.legend(handles = [l1,l2,l3])

    ax4=fig.add_subplot(gs[3])
    plt.plot(stableProfile, plume.lvl,label='vertical PM profile')
    ax4.set(xlabel='concentration [ug/kg]',ylabel='height [m]')
    ax4.fill_betweenx(plume.lvl, pmQ1, pmQ3, alpha=0.35,label='IQR')
    plt.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # plt.show()
    plt.savefig(plume.figdir + 'downwindAve/%s.pdf' %Case)

    plt.close()
    print('.....saved in: %s' %(plume.figdir + 'downwindAve/%s.pdf' %Case))

    #save some data for non-averaged comparison of fireline lengths
    if Case in plume.fireline_runs:
        Uprime = (csdict['u'][-1,:,:].T - U0).T
        FirelineUprime400m.append(Uprime[10,:])
        FirelineProfiles.append(stableProfile)

    #making a smoke plot for ENV Interagency presentation
    if Case == 'W5F6R1':
        fig = plt.figure(figsize=(14.5,4))
        ax1 = plt.gca()
        # ---u contours and colorbar
        im = ax1.imshow(csdict['pm25'][-1,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],vmax = np.max(csdict['pm25'][-1,:,:])/20,cmap=plt.cm.cubehelix_r)
        cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
        cbari.set_label('smoke concentration [ug/kg]')
        # ---non-filled vapor contours and colorbar
        # cntr = ax1.contour(PMcontours[-1,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
        ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
        ax1.axhline(y = zi[nCase], ls=':', c='darkgrey', label='BL height at ignition')
        ax1.set(ylabel='height AGL [m]')
        ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
        ax1.legend()
        plt.tight_layout()
        plt.show()
        plt.close()


# #Do dimenional analysis
# Pi1 = Phi * (g**0.5) / (Omega * (zi**0.5))
# Pi2 = zCL/zi
# fig = plt.figure(figsize=(12,6))
# plt.suptitle('DIMENSIONLESS ANALYSIS')
# plt.subplot(1,2,1)
# plt.suptitle('Colored by profile wind')
# plt.scatter(Pi1,Pi2, c=plume.read_tag('W',RunList),cmap=plt.cm.jet)
# # plt.gca().set(xlabel = 'zCL/r', ylabel='Wf*/Ua')
# plt.colorbar()
# plt.subplot(1,2,2)
# plt.suptitle('Colored by profile number')
# plt.scatter(Pi1,Pi2, c=plume.read_tag('R',RunList),cmap=plt.cm.tab20b)
# ax = plt.gca()
# for i, txt in enumerate(RunList):
#     ax.annotate(txt, (Pi1[i], Pi2[i]),fontsize=6)
# # plt.gca().set(xlabel = 'zCL/r', ylabel='Wf*/Ua')
# plt.colorbar()
# plt.subplots_adjust(top=0.85)
# plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.show()
# # plt.savefig(plume.figdir + 'DimAnalysis.pdf' )
# plt.close()

#now plot zCl as a function of w*
wStar = (g*Phi*zi/(Omega))**(1/3.)
slope, intercept, r_value, p_value, std_err = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])
print('Sum of residuals: %0.2f' %r_value)

fig = plt.figure(figsize=(12,6))
plt.suptitle('zCL=FCN(W*): R = %.2f' %r_value)
plt.subplot(1,2,1)
ax1=plt.gca()
plt.scatter(wStar, zCL,  c=plume.read_tag('W',RunList), cmap=plt.cm.jet)
plt.plot(wStar, intercept + slope*wStar, c='grey')
plt.colorbar(label='Ua wind speed [m/s]')
ax1.set(xlabel='w* [m/s]',ylabel='zCL [m]')
plt.subplot(1,2,2)
ax2=plt.gca()
plt.scatter(wStar, zCL,  c=width, cmap=plt.cm.RdYlGn_r)
plt.colorbar(label='total 2D burn intensity')
for i, txt in enumerate(RunList):
    ax2.annotate(txt, (wStar[i], zCL[i]),fontsize=6)
ax2.set(xlabel='w* [m/s]',ylabel='zCL [m]')
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'zCl_wStar.pdf' )
plt.show()


#Fireline length ANALYSIS
#------------Velocity Enhancement-----------
firelineSmoke = np.mean(FirelineProfiles)
plt.figure()
haxis = np.arange(dimX)*plume.dx
ax = plt.gca()
plt.title('FIRE-INDUCED WIND DYNAMCIS at 400 m AGL')
plt.plot(haxis, FirelineUprime400m[0], label='1 km')
plt.plot(haxis, FirelineUprime400m[1], label='2 km')
plt.plot(haxis, FirelineUprime400m[2], label='4 km')
# plt.axhline(y=0,xmin=0,xmax=16000,color='grey')
plt.xlabel('horizontal distance [m]')
plt.ylabel('Uprime [m/s]' )
plt.xlim([0,max(haxis)])
plt.tight_layout()
plt.legend()
# plt.savefig(plume.figdir + 'fireline/Uprime400m.pdf')
plt.show()
plt.close()

#------------Profile Comparison-----------
plt.figure(figsize=(12,4))
plt.title('DONWDIND CONCENTRATIONS')
plt.plot(FirelineProfiles[0]/np.max(FirelineProfiles[0]),plume.lvl)
plt.plot(FirelineProfiles[1]/np.max(FirelineProfiles[1]),plume.lvl)
plt.plot(FirelineProfiles[2]/np.max(FirelineProfiles[2]),plume.lvl)
plt.gca().set(xlabel='normalized concentration', ylabel='height [m]')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.legend()
# plt.savefig(plume.figdir + 'fireline/DownwindSmokeProfiles.pdf')
plt.show()
