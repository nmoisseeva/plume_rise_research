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
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import ma
from matplotlib import gridspec



#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

preIgnT = 1 		#boolean: whether to use pre-ignition temp profile or averaged upwind profile

#=================end of input===============

print('TRACER-BASED ANALYSIS OF LES PLUME DATA')
print('===================================')


param_dict = {'fire':[]}
profile_dict = {'wmin':[],'pm_wmin':[],'tilt':[]}
profile_dict['meta'] = 'wmin: profiles of minimum velocity; \
    pm_wmin: tracer profile at downwind location of minimum vertical velocity \
    tilt: grid number of wmin profile locations'

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
runCnt = len(RunList)

plume_inj = np.empty(runCnt) * np.nan
cumT = np.empty(runCnt)* np.nan
inv_cumT = np.empty(runCnt) * np.nan
plumeTilt = np.empty(runCnt)* np.nan
fireWidth = np.empty((runCnt))* np.nan
kinI = np.empty((runCnt))* np.nan
parcelHeat = np.empty(runCnt)* np.nan
T0 = np.empty((runCnt,len(plume.lvl))) * np.nan
U0 = np.empty((runCnt,len(plume.lvl))) * np.nan
PMprofiles = np.empty((runCnt,len(plume.lvl)))* np.nan
firelineTest = []

BLdict = {'Ua':np.empty(runCnt) * np.nan, \
            'Ti':np.empty(runCnt) * np.nan,\
            'zi': np.empty(runCnt)* np.nan,\
            'inversion': np.empty((runCnt))* np.nan}

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


    #check if this is a fireline length test case:
    if Case in plume.fireline_runs:
        firelineTest.append(nCase)

    #mask plume as being at last 50ppm---------------------------------
    pm = ma.masked_where(avedict['pm25'] <= 30, avedict['pm25'] )
    w = ma.masked_where(avedict['pm25'] <= 30, avedict['w'] )
    # pm = ma.masked_where(avedict['pm25'] <= avedict['pm25'].max()*0.001, avedict['pm25'] )
    # w = ma.masked_where(avedict['pm25'] <= avedict['pm25'].max()*0.001, avedict['w'] )

    PMmaxVidx = pm.argmax(0)
    xmax,ymax = np.nanargmax(PMmaxVidx), np.nanmax(PMmaxVidx)
    centerline = ma.masked_where(plume.lvl[PMmaxVidx] == 0, plume.lvl[PMmaxVidx])
    tilt = ymax/xmax                #tilt is calculated based on max centerline height, not injection height!

    if preIgnT:
        profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
        T0[nCase] = np.load(profpath)
    else:
        T0[nCase] = avedict['temp'][:,0]

    U0[nCase] = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    #charatertistics of plume temperature anomalyies---------------------
    diffT = ma.masked_where(avedict['pm25'] <= 30, (avedict['temp'].T-T0[nCase]).T)           #get temperature anomaly
    Tctr = np.array([diffT[i,ni] for ni, i in enumerate(PMmaxVidx)])
    Wctr = np.array([w[i,ni] for ni, i in enumerate(PMmaxVidx)])
    # Uctr = np.array([avedict['u'][i,ni] for ni, i in enumerate(PMmaxVidx)])



    #based slice location on temperature anomaly

    intXnan = np.argwhere(np.diff(np.sign(Tctr))).flatten()             #find where the interesepts are
    intX = intXnan[np.isfinite(Tctr[intXnan])]                          #get only the values in the plume
    intX = intX[intX > plume.wi]                                        #get rid of everything ahead of fireline
    intX = intX[PMmaxVidx[intX] > 0]                                    #remove any near-surface kinks
    if len(intX) < 2:
        sliceX = np.nanargmin(Wctr[plume.wi:]) + plume.wi
    else:
        sliceX = intX[1]
    sliceZ = PMmaxVidx[sliceX]

    #BL characteristics -------------------------------------------------


    dT = T0[nCase][1:]-T0[nCase][0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 3
    BLdict['Ti'][nCase] = T0[nCase][si+1]              #characteristic BL temperature
    zi_idx = np.argmax(gradT[si:]) + si                 #vertical level index of BL top
    BLdict['inversion'][nCase] = np.mean(dT[zi_idx:zi_idx+20])/plume.dz
    BLdict['zi'][nCase] = plume.dz * zi_idx
    BLdict['Ua'][nCase] = np.mean(U0[nCase][si:10])
    print('Wind speed near ground: %.1f m/s' %(BLdict['Ua'][nCase]))

    Uf = np.array(avedict['u'][int(zi_idx/2.),:] - U0[nCase][int(zi_idx/2.)])

    #fire behavior ------------------------------------------------------

    # #calculate total flux from ground and fire
    ignited = np.array([i for i in avedict['ghfx'] if i > 2])

    #convert to kinematic heat flux (from kW/m2)
    I = ignited * 1000/(1.2 * 1005)
    kinI[nCase] = sum(I)
    parcelHeat[nCase] =kinI[nCase] * plume.dx / BLdict['Ua'][nCase]

    #store data
    plumeTilt[nCase] = tilt
    fireWidth[nCase] = len(ignited)
    plume_inj[nCase] = plume.lvl[sliceZ]
    print('Plume top: %d m' %plume_inj[nCase])


    #estimate atmospheric heating and cumT --------------------------------
    #get cumulative T based on  delT
    delTplume = T0[nCase][int(si+1):sliceZ]-T0[nCase][si:sliceZ-1]
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

    dimZ, dimX = np.shape(avedict['w'])
    haxis = np.arange(dimX)*plume.dx
    pmLevels = np.geomspace(30,int(np.max(avedict['pm25']))/10,num=10)

    #plot contours
    PMcontours = ma.masked_where(avedict['pm25'] <= 30,avedict['pm25'] )
    fig = plt.figure(figsize=(18,10))
    plt.suptitle('%s' %Case)
    gs = fig.add_gridspec(ncols=3, nrows=2,height_ratios=[3,1])

    ax1=fig.add_subplot(gs[0,0])
    # plt.subplot(2,3,1,gs[0])
    plt.title('Time-averaged W')
    # ---w contours and colorbar
    im = ax1.imshow(avedict['w'], origin='lower',extent=[0,haxis[-1],0,plume.lvl[-1]], cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
    ax1.set_aspect('auto')
    cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    cbari.set_label('w $[m s^{-2}]$')
    ax1.set_xlabel('horizontal distance [m]')
    ax1.set_ylabel('height AGL [m]')
    ax1.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    ax1.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax1.axhline(y=plume_inj[nCase],ls='--', c='black',label='derived plume top')
    ax1.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax1.axvline(x = xmax*plume.dx, ls=':',c='darkgrey',label='tilt defintion')
    ax1.legend()
    # ---non-filled pm contours and colorbar
    cntr = ax1.contour(PMcontours, extent=[0,haxis[-1],0,plume.lvl[-1]],locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
    ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    # ---heat flux
    axh1 = ax1.twinx()
    axh1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh1.set_ylim([0,150])
    axh1.tick_params(axis='y', colors='red')
    ln = axh1.plot(haxis, avedict['ghfx'], 'r-')
    axh1.set_xlim([0,haxis[-1]])


    # plt.subplot(2,3,2,gs[0])
    ax2=fig.add_subplot(gs[0,1])
    plt.title('Time-averaged U')
    # ax2 = plt.gca()
    # ---u contours and colorbar
    # im = ax2.imshow((avedict['u'].T-U0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.RdBu_r, vmin=BLdict['Ua'][nCase]-5, vmax=BLdict['Ua'][nCase]+5)
    im = ax2.imshow((avedict['u'].T-U0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.RdBu_r, vmin=-4, vmax=4)
    ax2.set_aspect('auto')
    cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    cbari.set_label('u $[m s^{-2}]$')
    ax2.set_xlabel('horizontal distance [m]')
    ax2.set_ylabel('height AGL [m]')
    ax2.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    ax2.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax2.axhline(y=plume_inj[nCase],ls='--', c='black',label='plume injection height')
    ax2.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax2.axvline(x = xmax*plume.dx, ls=':',c='darkgrey',label='tilt defintion')
    # ---non-filled vapor contours and colorbar
    cntr = ax2.contour(PMcontours, extent=[0,dimX*plume.dx,0,plume.lvl[-1]], locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
    ax2.plot(np.arange(dimX)*plume.dx,centerline,ls='--', c='darkgrey' )
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

    ax3=fig.add_subplot(gs[0,2])
    plt.title('Time-averaged T')
    # ax3 = plt.gca()
    # ---u contours and colorbar
    im = ax3.imshow((avedict['temp'].T-T0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.cubehelix_r,vmin=-2, vmax=16)
    ax3.set_aspect('auto')
    cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    cbari.set_label('T [K]')
    ax3.set_xlabel('horizontal distance [m]')
    ax3.set_ylabel('height AGL [m]')
    ax3.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    ax3.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax3.axhline(y=plume_inj[nCase],ls='--', c='black',label='derived plume top')
    ax3.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax3.axvline(x = xmax*plume.dx, ls=':',c='darkgrey',label='tilt defintion')
    # ---non-filled vapor contours and colorbar
    cntr = ax3.contour(PMcontours, extent=[0,dimX*plume.dx,0,plume.lvl[-1]], locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
    ax3.plot(np.arange(dimX)*plume.dx,centerline,ls='--', c='darkgrey' )
    # ains = inset_axes(plt.gca(), width='40%', height='2%', loc=1)
    # cbar = fig.colorbar(cntr, cax=ains, orientation='horizontal')
    # cbar.set_label('PM2.5 mixing ratio $[ug/kg]$',size=8)
    # cbar.ax.tick_params(labelsize=8)
    # ---heat flux
    axh3 = ax3.twinx()
    ln = axh3.plot(np.arange(dimX)*plume.dx, avedict['ghfx'], 'r-')
    axh3.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh3.set_ylim([0,150])
    axh3.tick_params(axis='y', colors='red')
    axh3.set_xlim([0,dimX*plume.dx])


    ax4 = fig.add_subplot(gs[1,0])
    ax4.plot(haxis, Wctr)
    ax4.axhline(y=0, color='black')
    ax4.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax4.set_xlim([0,haxis[-1]])
    plt.xlabel('horizontal distance [m]')
    plt.ylabel('$w_{ctr}$ [m/s]')

    ax5 = fig.add_subplot(gs[1,1])
    ax5.plot(haxis, Uf)
    ax5.axhline(y=0, color='black')
    ax5.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax5.set_xlim([0,haxis[-1]])
    plt.xlabel('horizontal distance [m]')
    plt.ylabel('$u_{fire}$ [m/s]')

    ax6 = fig.add_subplot(gs[1,2])
    ax6.plot(haxis, Tctr)
    ax6.axhline(y=0, color='black')
    ax6.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    ax6.set_xlim([0,haxis[-1]])
    plt.xlabel('horizontal distance [m]')
    plt.ylabel('$T_{ctr}$ [K]')

    plt.subplots_adjust(top=0.85)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(plume.figdir + 'injectionH/ave%s.pdf' %Case)
    print('.....-->saved in: %s' %(plume.figdir + 'injectionH/ave%s.pdf' %Case))
    plt.close()


    #create actual aspect ratio plots---------------------------------------------
    fig = plt.figure(figsize=(14.5,12))

    plt.suptitle('%s' %Case)

    plt.subplot(311)
    ax1 = plt.gca()
    plt.title('Time-averaged W')
    im = ax1.imshow(avedict['w'], origin='lower',extent=[0,haxis[-1],0,plume.lvl[-1]], cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
    cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    cbari.set_label('w $[m s^{-2}]$')
    ax1.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',xlabel='horizontal distance [m]',ylabel='height AGL [m]')
    ax1.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    ax1.legend()
    cntr = ax1.contour(PMcontours, extent=[0,haxis[-1],0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.8)
    axh1 = ax1.twinx()
    axh1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh1.tick_params(axis='y', colors='red')
    axh1.set(ylim=[0,150],xlim=[0,haxis[-1]])
    ln = axh1.plot(haxis, avedict['ghfx'], 'r-')

    plt.subplot(312)
    ax2 = plt.gca()
    plt.title('Time-averaged U anomaly')
    im = ax2.imshow((avedict['u'].T-U0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.RdBu_r, vmin=-4, vmax=4)
    cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    cbari.set_label('u $[m s^{-2}]$')
    ax2.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',xlabel='horizontal distance [m]',ylabel='height AGL [m]')
    ax2.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax2.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )

    ax2.legend()
    cntr = ax2.contour(PMcontours, extent=[0,haxis[-1],0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.8)
    axh2 = ax2.twinx()
    axh2.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh2.tick_params(axis='y', colors='red')
    axh2.set(ylim=[0,150],xlim=[0,haxis[-1]])
    ln = axh2.plot(haxis, avedict['ghfx'], 'r-')

    plt.subplot(313)
    ax3 = plt.gca()
    plt.title('Time-averaged T anomaly')
    ax3 = plt.gca()
    im = ax3.imshow((avedict['temp'].T-T0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.cubehelix_r,vmin=-2, vmax=16)
    cbari =fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    cbari.set_label('w $[m s^{-2}]$')
    ax3.set(xlim=[0,haxis[-1]],ylim=[0,plume.lvl[-1]],aspect='equal',xlabel='horizontal distance [m]',ylabel='height AGL [m]')
    ax3.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax3.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    ax3.legend()
    cntr = ax3.contour(PMcontours, extent=[0,haxis[-1],0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.8)
    axh3 = ax3.twinx()
    axh3.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh3.tick_params(axis='y', colors='red')
    axh3.set(ylim=[0,150],xlim=[0,haxis[-1]])
    ln = axh3.plot(haxis, avedict['ghfx'], 'r-')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(plume.figdir + 'fixedAspectAverages/ave%s.pdf' %Case)
    print('.....-->saved in: %s' %(plume.figdir + 'fixedAspectAverages/ave%s.pdf' %Case))
    # plt.show()
    plt.close()


#---------------------calculations over all plumes--------------
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)
Ftag = np.array([i for i in plume.read_tag('F',RunList)])
#WHAT LEVEL TO PICK FOR defining characteristic U
# print('Variability of normalized intensity for given U level:')
# print('%d.2' %np.std(parcelHeat))

#-------------subplots of tilt predictors-----------
#create scatter plots of tilts vs other variables
plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
plt.title('FIRELINE INTENSITY vs TILT')
ax1 = plt.scatter(parcelHeat,plumeTilt,c=BLdict['Ua'],cmap=plt.cm.viridis)
plt.xlabel('parcel heat [K m]')
plt.ylabel('plume tilt')
plt.colorbar(ax1,label='windspeed [m/s]')

plt.subplot(1,3,2)
plt.title('WIND vs TILT')
ax2 = plt.scatter(BLdict['Ua'],plumeTilt/(kinI/fireWidth),c=BLdict['zi'],cmap=plt.cm.plasma)
plt.scatter(BLdict['Ua'][firelineTest], plumeTilt[firelineTest]/(kinI[firelineTest]/fireWidth[firelineTest]), c='green',marker='^')
# ax2 = plt.scatter(BLdict['Ua'],plumeTilt/parcelHeat,c=BLdict['zi'],cmap=plt.cm.plasma)

plt.xlabel('near-surface wind (m/s)')
plt.ylabel('plume tilt/(kinI/fireWidth) [$m_x \cdot s \cdot K^{-1} \cdot m_z^{-1}$]')
plt.colorbar(ax2, label='BL height [m]')

plt.subplot(1,3,3)
plt.title('WIDTH vs TILT')
ax3 = plt.scatter(fireWidth,plumeTilt,c=BLdict['Ua'],cmap=plt.cm.viridis)
plt.xlabel('fireline width [#grids]')
plt.ylabel('plume tilt')
plt.colorbar(ax3,label='windspeed [m/s]')

plt.tight_layout()
plt.savefig(plume.figdir + 'tilt_subplots.pdf')
# plt.show()
# plt.close()


#attempt at characteristic vertical velocity
c = BLdict['Ua']**(-2.) * BLdict['inversion']**(-1.) * BLdict['zi']**(-1.)
fromplot=plumeTilt * BLdict['Ua']
fromtest = c * BLdict['Ua']**(2.) * kinI/fireWidth
plt.scatter(fromplot,fromtest,c=kinI/fireWidth)
regV = np.poly1d(np.polyfit(fromplot,fromtest,1))
plt.xlim([0,10])
plt.ylim([0,10])
print(regV)
#-------------subplots for eatch R run-----------
plt.figure(figsize=(18,12))
for R in set(Rtag):
    plt.subplot(2,3,R+1)
    plt.title('R%s RUNS' %R, color='C%s' %R)
    ax1 = plt.scatter(parcelHeat[Rtag==R],cumT[Rtag==R],marker='o', c=plume.read_tag('W',RunList)[Rtag==R], cmap=plt.cm.viridis)
    for i, txt in enumerate(plume.read_tag('F',RunList)[Rtag==R]):
        plt.gca().annotate(txt, (parcelHeat[Rtag==R][i]+20,cumT[Rtag==R][i]+20), fontsize=9)
    plt.xlabel('normalized fireline intensity [$K m$]')
    plt.ylabel('cumulative temperature [K m]')
    plt.colorbar(ax1,label='windspeed [m/s]')
    plt.xlim([0,7000])
    plt.ylim([0,400])
plt.tight_layout()
plt.savefig(plume.figdir + 'Rn_sublots.pdf')
plt.show()
plt.close()

#-------------main plot with CumT-----------
plt.figure()
plt.title('NORMALIZED FIRELINE INTENSITY vs CumT')
ax = plt.gca()
for R in set(Rtag):
    regR = np.poly1d(np.polyfit(parcelHeat[Rtag==R],cumT[Rtag==R],1))
    print('R%s profiles: regression fit equation: %s ' %(R,regR))
    ax.scatter(parcelHeat[Rtag==R],cumT[Rtag==R], c='C%s'%R ,label='R%s' %R)
    plt.plot(parcelHeat[Rtag==R], regR(parcelHeat[Rtag==R]), 'C%s' %R)
ax.scatter(parcelHeat[firelineTest], cumT[firelineTest], c='green',marker='^')
plt.xlabel('normalized fireline intensity [$K m$]')
plt.ylabel('cumulative temperature [K m]')
plt.legend()
plt.tight_layout()
plt.savefig(plume.figdir + 'parcelHeat_cumT.pdf')
plt.show()
plt.close()

#-------------plot initial profiles---------------
plt.figure()
plt.title('INITIAL ATMOSPHERIC PROFILES')
leg_handles = []
for R in set(Rtag):
    for Case in T0[Rtag==R]:
        lR = plt.plot(Case, plume.lvl, color='C%s' %R, label='R%s' %R)
    leg_handles.extend(lR)
plt.xlabel('potential temperature [K]')
plt.ylabel('height [m]')
plt.legend(handles=leg_handles)
plt.savefig(plume.figdir + 'T0profiles.pdf')
plt.show()


#how to estaimte bounary curvature
# At each boundary point, we calculate the boundary curvature by fitting a circle to that boundary point and the two points that are 10 boundary points away from it. The magnitude of the boundary curvature is then defined as the reciprocal of the radius of that circle. If the midpoint of the two points 10 boundary points away is inside the cell, the curvature is defined as positive, otherwise it is defined as negative. For visualization, the curvature is smoothed over 3 boundary points and 3 frames, and the color scale is cut off at a maximum curvature magnitude.
#tracting boundaries: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3305346/

# plt.figure(figsize=(12,6))
# clr = plt.cm.plasma(plt.Normalize()(parcelHeat))
# # clr[..., -1] = plt.Normalize()(fireLine[:,0])
# for nCase,Case in enumerate(RunList):
#     plt.subplot(1,2,2)
#     plt.title('NORMALIZED VERTICAL Q/Qmax PROFILES')
#     plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_inj[nCase],c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile ')
#     # plt.colorbar(label='fireline intensity ')
#     plt.ylim([0,2])
#     plt.subplot(1,2,1)
#     plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
#     plt.plot(Qprofiles[nCase,:]/Ua[nCase],plume.lvl/plume_inj[nCase], c=clr[nCase] )
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
#     plt.plot(Qprofiles[nCase,:]/(np.max(Qprofiles[nCase,:])),plume.lvl/plume_inj[nCase],c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile ')
#     # plt.colorbar(label='fireline intensity ')
#     plt.ylim([0,2])
#     plt.subplot(1,2,1)
#     plt.title('NORMALIZED VERTICAL Q/Umean PROFILES')
#     plt.plot(Qprofiles[nCase,:]/Ua[nCase],plume.lvl/plume_inj[nCase], c=clr[nCase] )
#     plt.ylabel('normalized height')
#     plt.xlabel('normalized $H_2O$ profile [g s / kg m]')
#     plt.ylim([0,2])
# plt.savefig(plume.figdir + 'normProfsHFX.pdf')
# # plt.show()
# plt.close()
