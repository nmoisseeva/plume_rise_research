# January 2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter



#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time


#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
# RunList = ['W4F7R1']
runCnt = len(RunList)


#take a vertical cross-section (ensure it's stable)
#calculate CDF of cumT (delT*dz)
#calculate temperature gradient -
#map temperature gradient to CDF


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
    zi = plume.get_zi(T0)

    ctrZidx = pm.argmax(0)
    ctrXidx = pm.argmax(1)
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])
    for nZ in range(dimZ):
        if pmCtr[ctrXidx[nZ]] < pm[nZ,ctrXidx[nZ]]:
            pmCtr[ctrXidx[nZ]] = pm[nZ,ctrXidx[nZ]]


    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
    centerline = ma.masked_where(plume.lvl[ctrZidx] == 0, plume.lvl[ctrZidx])
    tilt = ymax/xmax                #tilt is calculated based on max centerline height


    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 81, 2) # window size 51, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]

    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.mean(stablePM,1)
    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)



    #PLOTTING =========================================================
    haxis = np.arange(dimX)*plume.dx
    maxPM = int(np.max(csdict['pm25']))
    pmLevels = np.geomspace(30,maxPM/10,num=10)
    PMcontours = ma.masked_where(csdict['pm25'] <= 30,csdict['pm25'])

    fig = plt.figure(figsize=(22,8))
    gs = fig.add_gridspec(ncols=2, nrows=2,width_ratios=[6,1])
    print('.....creating vertical crossection of U + PM2.5')
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
    ax1.axhline(y = zi, ls=':', c='darkgrey', label='BL height at ignition')
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
    plt.plot(haxis, pmCtr, label='concentration along centerline', color='C1')
    ax3.fill_between(haxis[1:], 0, 1, where=stablePMmask, color='grey', alpha=0.4, transform=ax3.get_xaxis_transform(), label='averaging window')
    ax3.set(xlim=[0,dimX*plume.dx],xlabel='horizontal distance [m]',ylabel='concentration [ug/kg]' )
    plt.legend()

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
