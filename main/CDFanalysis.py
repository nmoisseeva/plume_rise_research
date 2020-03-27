# February 2020

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
# RunList = ['W4F7R4L4']

runCnt = len(RunList)
g = 9.81
si = 3

#save variables for dimensional analysis
zi = np.empty((runCnt)) * np.nan            #BL height
zCL = np.empty((runCnt)) * np.nan           #injection height
Omega = np.empty((runCnt)) * np.nan         #cumulative BL resistance
Ir = np.empty((runCnt)) * np.nan            #cumulative fire intensity
xMax = np.empty((runCnt)) * np.nan          #horiztonal location of certerline peak
Ua = np.empty((runCnt)) * np.nan                #ambient wind


for nCase,Case in enumerate(RunList):
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)
    csdict = plume.load_CS_prep_Profiles(Case)
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    #define heat source------------------------
    masked_flux = ma.masked_less_equal(csdict['ghfx2D'], 0)
    cs_flux = np.nanmean(masked_flux,1)                         #get cross section for each timestep
    fire = []
    fmax = np.argmax(cs_flux,axis=1)                            #get max for each timestep
    for nP, pt in enumerate(fmax[plume.ign_over:]):             #excludes steps containing ignition
        subset = cs_flux[plume.ign_over+nP,pt-plume.wi:pt+plume.wf]     #set averaging window around a maximum
        fire.append(subset)

    meanFire = np.nanmean(fire,0)
    ignited = np.array([i for i in meanFire if i > 0.5]) * 1000 / ( 1.2 * 1005)
    Ir[nCase] = np.trapz(ignited, dx = plume.dx)

    #mask plume with cutoff value---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= plume.PMcutoff, csdict['pm25'][-1,:,:] )
    zi[nCase] = plume.get_zi(T0)
    si = 3

    #locate centerline
    ctrZidx = pm.argmax(0)
    ctrZidx[:fmax[-1]] = 0
    ctrXidx = pm.argmax(1)
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])

    for nZ in range(dimZ):
        if pmCtr[ctrXidx[nZ]] < pm[nZ,ctrXidx[nZ]]:
            pmCtr[ctrXidx[nZ]] = pm[nZ,ctrXidx[nZ]]

    #define downwind distribution
    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
    if fmax[-1] > xmax:
        print('\033[93m' + 'Fire ahead of overshoot!!!' +  '\033[0m')
    centerline = ma.masked_where(plume.lvl[ctrZidx] == 0, plume.lvl[ctrZidx])
    smoothCenterline = savgol_filter(centerline, 31, 3) # window size 101, polynomial order 3
    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 101, 3) # window size 101, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.05 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]
    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.mean(stablePM,1)
    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])
    zCLidx = int(np.mean(ctrZidx[1:][stablePMmask]))

    xMax[nCase] = ymax/xmax
    #get fintal profile weights by calculating proportional area for each level
    subArea = stableProfile * plume.dx
    totalArea = sum(subArea.data)
    weights = subArea/totalArea



    #dimensional analysis variables ---------------------------
    dT = T0[1:]-T0[0:-1]
    Ua[nCase] = np.mean(U0[si:zCLidx])
    # cumT = np.trapz(np.cumsum(dT[si+1:])* weights[si+1:-1], dx=plume.dz)
    cumT = np.trapz(dT[si+1:zCLidx], dx = plume.dz)       #this is what we typically use

    # Omega[nCase] = np.nansum(cumT)
    Omega[nCase] = cumT


    #
    #
    # #PLOTTING =========================================================
    # haxis = np.arange(dimX)*plume.dx
    # maxPM = int(np.max(csdict['pm25']))
    # pmLevels = np.geomspace(plume.PMcutoff,maxPM/10,num=10)
    # PMcontours = ma.masked_where(csdict['pm25'] <= plume.PMcutoff,csdict['pm25'])
    #
    # fig = plt.figure(figsize=(22,8))
    # gs = fig.add_gridspec(ncols=2, nrows=2,width_ratios=[6,1])
    # plt.suptitle('%s' %Case)
    #
    # ax1=fig.add_subplot(gs[0])
    # axh1=ax1.twinx()
    # # ---u contours and colorbar
    # im = ax1.imshow((csdict['u'][-1,:,:].T - U0).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)
    # cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    # cbari.set_label('ralative horizontal velocity $[m s^{-1}]$')
    # # ---non-filled vapor contours and colorbar
    # cntr = ax1.contour(PMcontours[-1,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
    # ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    # ax1.axhline(y = zi[nCase], ls=':', c='darkgrey', label='BL height at ignition')
    # ax1.set(ylabel='height AGL [m]')
    # ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
    # ax1.legend()
    # # ---heat flux
    # ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][-1,:], 'r-')
    # axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
    # axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])
    # axh1.tick_params(axis='y', colors='red')
    #
    #
    # ax2=fig.add_subplot(gs[1])
    # fim = ax2.imshow(csdict['ghfx2D'][-1,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
    # cbarif = fig.colorbar(fim, orientation='horizontal')
    # cbarif.set_label('heat flux [$kW / m^2$]')
    # ax2.set(xlabel='x distance [m]',ylabel='y distance [m]',aspect='equal')
    #
    #
    # ax3=fig.add_subplot(gs[2])
    # l1, = plt.plot(haxis, pmCtr, label='concentration along centerline', color='C1')
    # l2 = ax3.fill_between(haxis[1:], 0, 1, where=stablePMmask, color='grey', alpha=0.4, transform=ax3.get_xaxis_transform(), label='averaging window')
    # ax3.set(xlim=[0,dimX*plume.dx],xlabel='horizontal distance [m]',ylabel='concentration [ug/kg]' )
    # ax32 = ax3.twinx()
    # l3, = plt.plot(haxis,smoothCenterline, label='smoothed centerline height ', color='C2',linestyle=':')
    # ax32.set(xlim=[0,dimX*plume.dx],ylim=[0,2000], ylabel='height [m]' )
    # plt.legend(handles = [l1,l2,l3])
    #
    # ax4=fig.add_subplot(gs[3])
    # plt.plot(stableProfile, plume.lvl,label='vertical PM profile')
    # ax4.set(xlabel='concentration [ug/kg]',ylabel='height [m]')
    # ax4.fill_betweenx(plume.lvl, pmQ1, pmQ3, alpha=0.35,label='IQR')
    # plt.legend()
    #
    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    # # plt.show()
    # # plt.savefig(plume.figdir + 'downwindAve/%s.pdf' %Case)
    #
    # plt.close()
    # # print('.....saved in: %s' %(plume.figdir + 'downwindAve/%s.pdf' %Case))
    #


##now plot zCl as a function of w*
wStar = (g*Ir*zi/(Omega))**(1/3.)
Nf = np.sqrt(g*Ir/(Omega * zi * Ua))

# slope, intercept, r_value, p_value, std_err = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])
# print('Sum of residuals: %0.2f' %r_value)

fig = plt.figure(figsize=(12,6))
# plt.suptitle('zCL=FCN(W*): R = %.2f' %r_value)
plt.subplot(1,2,1)
ax1=plt.gca()
plt.scatter(wStar, zCL,  c=plume.read_tag('W',RunList), cmap=plt.cm.jet)
# plt.plot(wStar, intercept + slope*wStar, c='grey')
plt.colorbar(label='Ua wind speed [m/s]')
ax1.set(xlabel='$w*',ylabel='zCL')
plt.subplot(1,2,2)
ax2=plt.gca()
plt.scatter(wStar, zCL,  c=zi, cmap=plt.cm.RdYlGn_r)
plt.colorbar(label='zi')
# for i, txt in enumerate(RunList):
#     ax2.annotate(txt, (wStar[i], zCL[i]),fontsize=6)
# ax2.set(xlabel='w* [m/s]',ylabel='zCL [m]')
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig(plume.figdir + 'zCl_wStar.pdf' )
plt.show()
plt.close()


fig = plt.figure(figsize=(12,6))
# plt.suptitle('zCL=FCN(W*): R = %.2f' %r_value)
plt.subplot(1,2,1)
ax1=plt.gca()
plt.scatter(Ir, Omega,  c=plume.read_tag('W',RunList), cmap=plt.cm.jet)
# plt.plot(wStar, intercept + slope*wStar, c='grey')
plt.colorbar(label='Ua wind speed [m/s]')
ax1.set(xlabel='Ir',ylabel='Omega')
plt.subplot(1,2,2)
ax2=plt.gca()
plt.scatter(Ir, Omega,  c=zi, cmap=plt.cm.RdYlGn_r)
plt.colorbar(label='zi')
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig(plume.figdir + 'zCl_wStar.pdf' )
# plt.show()
plt.close()

#plot "fire frequency" vs tilt
slope, intercept, r_value, p_value, std_err = linregress(Nf[np.isfinite(Nf)],xMax[np.isfinite(Nf)])
# print('Sum of residuals: %0.2f' %r_value)
