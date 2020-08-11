# May 2020
#This code plots CWI smoke on last frame and applies filter to determine zCL_true

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.interpolate import interp1d




#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time
#=================end of input===============

RunList =   [i for i in plume.tag if i not in plume.exclude_runs]
# RunList = ['W5F2R6T','W5F3R6T','W5F4R6T','W5F5R6T','W5F6R6T','W5F7R6T','W5F10R6T','W5F12R6T','W5F13R6T', 'W3F7R6T','W4F7R6T','W6F7R6T','W7F7R6T','W8F7R6T','W9F7R6T','W10F7R6T','W11F7R6T','W12F7R6T']
# RunList = ['W5F7R0','W5F7R1','W5F7R2','W5F7R3','W5F7R4','W5F7R5T','W5F7R6T']

runCnt = len(RunList)
g = 9.81

#set up interpolated vertical profile with 5m vertical step
interpZ = np.arange(0, plume.lvltall[-1], plume.zstep)
si = int(plume.zs/plume.zstep)

#storage for variables
zi = np.empty((runCnt)) * np.nan                #BL height (m)
zCL = np.empty((runCnt)) * np.nan               #smoke injection height (m)
Phi = np.empty((runCnt)) * np.nan               #cumulative fire heat  (K m^2/s)
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature (Km) - the questionable denominator term
sounding = np.empty((runCnt,len(interpZ))) * np.nan         #storage for interpolated soundings
gradT0interp = np.empty((runCnt,len(interpZ)-1)) * np.nan   #storage for temperature gradient
zCLerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
profile = np.empty((runCnt,len(interpZ))) * np.nan
quartiles = np.empty((runCnt,len(interpZ),2)) * np.nan

#======================repeat main analysis for all runs first===================
#loop through all LES cases
for nCase,Case in enumerate(RunList):
    #exclude outlier runs that are undealt with
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)

    csdict = plume.load_CS_prep_Profiles(Case)      #load data cross-wind integrated smoke and all other data
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')    #load initial temperature profile
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')    #load intial wind profile

    #create an interpolated profile of temperature
    if Case[-1:]=='T' or Case[-1:]=='E':
        levels = plume.lvltall
    else:
        levels=plume.lvl
    interpT= interp1d(levels,T0,fill_value='extrapolate')
    T0interp = interpT(interpZ)

    #mask plume with cutoff value---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])     #get shape of data
    zi[nCase] = plume.get_zi(T0)                    #calculate BL height
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= plume.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells

    #locate centerline
    ctrZidx = pm.argmax(0)                          #locate maxima along height
    ctrXidx = pm.argmax(1)                          #locate maxima downwind
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])    #get concentration along the centerline


    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)           #get location of maximum centerline height
    centerline = ma.masked_where(plume.lvltall[ctrZidx] == 0, plume.lvltall[ctrZidx])               #make sure centerline is only calculated inside the plume
    centerline.mask[:int(1000/plume.dx)] = True
    # smoothCenterline = savgol_filter(centerline, 51, 3)             # smooth centerline height (window size 31, polynomial order 3)

    filter_window = max(int(plume.read_tag('W',[Case])*10+1),51)
    smoothCenterline = savgol_filter(centerline, filter_window, 3)             # smooth centerline height (window size 31, polynomial order 3)


    #calculate concentration changes along the centerline
    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    # smoothPM = savgol_filter(dPMdX, 101, 3) # window size 101, polynomial order 3
    smoothPM = savgol_filter(dPMdX, filter_window, 3) # window size 101, polynomial order 3

    # stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.05 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                            abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                            nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                            nX > np.nanargmax(smoothPM) and\
                            nX > np.nanargmax(centerline) +10 and\
                            centerline[nX] < plume.lvltall[-1]-200 and \
                            nX > np.nanargmax(smoothCenterline)+10 else \
                            False for nX in range(dimX-1) ]
    if sum(stablePMmask) == 0:
        stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                                nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                nX > np.nanargmax(centerline) +10 and\
                                nX > np.nanargmax(smoothPM) else\
                                # nX > np.nanargmax(smoothCenterline) else \
                                False for nX in range(dimX-1) ]
    if sum(stablePMmask) == 0:
        stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                                nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                nX > np.nanargmax(smoothPM) else\
                                # nX > np.nanargmax(smoothCenterline) else \
                                False for nX in range(dimX-1) ]
    if Case=='W5F4R6T':
        x = np.array(stablePMmask)
        x[:int(6000/plume.dx)] = False
        stablePMmask = list(x)
    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.median(stablePM,1)

    interpT= interp1d(levels,T0,fill_value='extrapolate')
    T0interp = interpT(interpZ)

    profile[nCase,:] = interp1d(levels,stableProfile,fill_value='extrapolate')(interpZ)

    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)
    quartiles[nCase,:,0] = interp1d(levels,pmQ1,fill_value='extrapolate')(interpZ)
    quartiles[nCase,:,1] = interp1d(levels,pmQ3,fill_value='extrapolate')(interpZ)

    #define heat source ------------------------
    masked_flux = ma.masked_less_equal(csdict['ghfx2D'],1)    #mask empty fire heat flux cells
    cs_flux = np.nanmean(masked_flux,1)                         #get mean cross section for each timestep
    fire = []                                                   #create storage arrage
    fxmax = np.argmax(cs_flux,axis=1)                           #get location of max heat for each timestep
    for nP, pt in enumerate(fxmax[plume.ign_over:]):            #excludes steps containing ignition
        subset = cs_flux[plume.ign_over+nP,pt-plume.wi:pt+plume.wf]     #set averaging window around a maximum
        fire.append(subset)

    meanFire = np.nanmean(fire,0)                               #calculate mean fire cross section
    ignited = np.array([i for i in meanFire if i > 0.5])        #consider only cells that have heat flux about 500 W/m2

    Phi[nCase] = np.trapz(ignited, dx = plume.dx) * 1000 / ( 1.2 * 1005)    #calculate Phi by integrating kinematic heat flux along x (Km2/s)

    #calculate injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(interpZ - zCL[nCase]))

    dT = (T0[1:]-T0[0:-1])/plume.dz                                        #calculate potential temperature change (K)

    sounding[nCase,:] = T0interp
    dTinterp = (T0interp[1:] - T0interp[0:-1])/plume.zstep
    gradT0interp[nCase,:] = dTinterp
    Omega[nCase] = np.trapz(dTinterp[si:zCLidx], dx = plume.zstep)     #calculate the denominator term by integrating temperature change with height, excluding surface layer (Km)



    #PLOTTING =========================================================
    cropX = int(dimX*0.75)
    axMax = cropX * plume.dx
    haxis = np.arange(cropX)*plume.dx
    PMppm = pm/1000.                                        #smoke concentration in ppm
    maxPM = int(np.max(PMppm))

    fig = plt.figure(figsize=(11,5))
    gs = fig.add_gridspec(ncols=2, nrows=2,width_ratios=[4,1])
    plt.suptitle('%s' %Case)

    ax1=fig.add_subplot(gs[0])
    axh1=ax1.twinx()
    # ---cwi smoke  and colorbar
    im = ax1.imshow(PMppm[:,:cropX], origin='lower', extent=[0,axMax,0,levels[-1]],cmap=plt.cm.cubehelix_r,vmin=0, vmax=maxPM/10)
    cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
    cbari.set_label('CWI smoke $[ppm]$')
    ax1.plot(haxis,centerline[:cropX],ls='--', c='dimgrey',label='plume centerline' )
    ax1.axhline(y = zi[nCase], ls=':', c='darkgrey', label='BL height at ignition')
    ax1.set(ylabel='height [m]')
    ax1.set(xlim=[0,axMax],ylim=[0,levels[-1]],aspect='equal')
    ax1.legend()
    # ---heat flux
    ln = axh1.plot(haxis, csdict['ghfx'][-1,:cropX], 'r-')
    axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
    axh1.set(xlim=[0,axMax],ylim=[0,150])
    axh1.tick_params(axis='y', colors='red')
    ax1.text(0.02, 0.9, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, weight='bold')


    ax2=fig.add_subplot(gs[1])
    fim = ax2.imshow(csdict['ghfx2D'][-1,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
    cbarif = fig.colorbar(fim, orientation='vertical')
    cbarif.set_label('heat flux [$kW / m^2$]')
    ax2.set(xlabel='x distance [m]',ylabel='y distance [m]',aspect='equal')
    ax2.text(0.1, 0.93, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, weight='bold')


    ax3=fig.add_subplot(gs[2])
    l1 = ax3.fill_between(haxis, 0, 1, where=stablePMmask[:cropX], color='grey', alpha=0.4, transform=ax3.get_xaxis_transform(), label='averaging window')
    ax3.set(xlim=[0,axMax],ylim=[0,3200], ylabel='height [m]',xlabel='distance [m]')
    l3, = ax3.plot(haxis,smoothCenterline[:cropX], label='smoothed centerline height ', color='C2')
    l2, = ax3.plot(haxis,centerline[:cropX], label='raw centerline height', color='C4',linestyle=':')
    ax32 = ax3.twinx()
    ax32.set(xlim=[0,axMax],xlabel='distance [m]')
    l4, = plt.plot(haxis, pmCtr[:cropX]/1000, label='concentration gradient', color='C1',linewidth=1)
    plt.legend(handles = [l1,l2,l3,l4])
    ax3.text(0.02, 0.93, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, weight='bold')
    ax32.tick_params(axis='y',colors='C1')
    ax32.set_ylabel('concentration gradient [ppm]', color='C1')

    ax4=fig.add_subplot(gs[3])
    plt.plot(stableProfile/1000, levels,label=' PM profile')
    ax4.set(xlabel='CWI concentration [ppm]',ylabel='height [m]')
    ax4.fill_betweenx(levels, pmQ1/1000, pmQ3/1000, alpha=0.35,label='IQR')
    ax4.axhline(y = zCL[nCase], ls='--', c='black', label='z$_{CL}$')
    ax4.text(0.1, 0.93, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, weight='bold')

    plt.legend()
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # plt.show()
    plt.savefig(plume.figdir + 'CWIzCL/zcl%s.pdf' %Case)
    plt.close()
    print('.....saved in: %s' %(plume.figdir + 'CWIzCL/zcl%s.pdf' %Case))


BLfrac = 0.75           #fraction of BL height to set zs at
mf, bf = 0.9243, 115.059
C = 1.0087
thetaS = np.empty((runCnt)) * np.nan               #smoke injection height (m)
thetaCL = np.empty((runCnt)) * np.nan               #smoke injection height (m)
OmegaOver = np.empty((runCnt)) * np.nan               #smoke injection height (m)
OmegaUnder = np.empty((runCnt)) * np.nan               #smoke injection height (m)
zCLP = np.empty((runCnt)) * np.nan
Gamma = np.empty((runCnt)) * np.nan
zCLGuess = np.empty((runCnt)) * np.nan
zMax = np.empty((runCnt)) * np.nan
zMaxGuess = np.empty((runCnt)) * np.nan

exclude = ['W5F4R6TE','W5F13R6TE','W5F12R5TE','W5F4R5TE','W5F1R1','W5F13R7T','W5F7R8T','W5F13R5TE']

for nCase,Case in enumerate(RunList):
    if Case in exclude:
        continue
    else:
        zS = zi[nCase]*BLfrac
        zsidx = np.argmin(abs(interpZ - zS))
        ziidx = np.argmin(abs(interpZ - zi[nCase]))
        zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
        zCLidxP = np.argmax(profile[nCase])
        thetaS[nCase] = sounding[nCase,zsidx]
        thetaCL[nCase] = sounding[nCase,zCLidx]

        #far zCL estimate
        Tau = 1/ np.sqrt(g*Omega[nCase]/(thetaS[nCase] * (zCL[nCase]-zS)))
        wStar = C*Tau*((g*Phi[nCase]*(zCL[nCase]-zS))/(thetaS[nCase]*zi[nCase]))**(1/3.)
        zCLGuess[nCase] = mf*(wStar+zS) + bf

        dPM = profile[nCase,1:] - profile[nCase,:-1]
        cutoff = dPM[zCLidxP+5:][abs(dPM[zCLidxP+5:])<np.max(dPM)*0.05][0]
        pmTopFind = np.where(dPM==cutoff)
        if len(pmTopFind)==0:
            pmTopidx = len(interpZ)
        else:
            pmTopidx = pmTopFind[0][0]
        OmegaUnder[nCase] = sounding[nCase,zCLidxP] - sounding[nCase,zsidx]
        OmegaOver[nCase] = sounding[nCase,pmTopidx] - sounding[nCase,zCLidxP]
        zMax[nCase] = interpZ[pmTopidx]
        zCLP[nCase] = interpZ[zCLidxP]

        #gamma fit
        # fittopidx = np.nanargmin(abs(interpZ - 2700))
        baselinedTheta = sounding[nCase,zCLidxP:zCLidxP+50] - sounding[nCase,zsidx]
        GammaFit = linregress(interpZ[zCLidxP:zCLidxP+50],baselinedTheta)
        Gamma[nCase] = GammaFit[0]

        #get top prediction
        overshoot = (zCLP[nCase] - zS)*Gamma[nCase]
        dTabove = sounding[nCase,zCLidxP+1:] - sounding[nCase,zCLidxP:-1]
        OmegaOverGuess = np.cumsum(dTabove)
        zMaxidx= np.argmin(abs(OmegaOverGuess - overshoot)) + zCLidxP
        zMaxGuess[nCase] = interpZ[zMaxidx]

        plt.figure()
        plt.title('%s' %Case)
        plt.plot(profile[nCase,:]/1000,interpZ,label=' PM median profile')
        ax = plt.gca()
        ax.set(xlabel='CWI concentration [ppm]',ylabel='height [m]')
        ax.fill_betweenx(interpZ, quartiles[nCase,:,0]/1000,quartiles[nCase,:,1]/1000, alpha=0.35,label='IQR')
        ax.axhline(y = interpZ[zCLidxP], ls='--', c='black', label='z$_{CL}$ profile')
        # ax.axhline(y = zCLGuess[nCase], ls='--', c='red', label='z$_{CL} LES$')
        ax.axhline(y = zMax, ls=':', c='purple',label='$z_{max}$ true' )
        ax.axhline(y = zMaxGuess, ls=':', c='red',label='$z_{max}$ true' )

        plt.legend()
        plt.savefig(plume.figdir + 'distribution/pmProf%s.pdf' %Case)
        plt.close()

errorMax = zMax-zMaxGuess
predictor = (zCLP-zi*BLfrac)*Gamma
topFit = linregress(predictor[np.isfinite(predictor)],OmegaOver[np.isfinite(OmegaOver)])

plt.figure(figsize=(10,5))
plt.subplot(121)
plt.title('PREDICTING DISTRIBUTION TOP: R=%.2f' %topFit[0])
plt.scatter((zCLP-zi*BLfrac)*Gamma,OmegaOver)
plt.gca().set(xlabel='$z^\prime \gamma$ [K]', ylabel='$\\theta_{top} - \\theta_{CL}$ [K]',aspect='equal',xlim=[0,20],ylim=[0,20])

plt.subplot(122)
plt.title(r'ERROR BOXPLOT: (true - model)')
plt.boxplot(errorMax[np.isfinite(errorMax)])
plt.tight_layout()
plt.savefig(plume.figdir + 'distribution/OmegaOver.pdf' )
plt.close()

# plt.xlim([0,0.0002])

# for i, txt in enumerate(RunList):
#     plt.gca().annotate(txt, (predictor[i], OmegaOver[i]),fontsize=6)

# np.savetxt('profiles.csv',sounding.T,fmt='%.2f',delimiter=',', header = str(plume.read_tag('R',RunList)))
# plt.figure(figsize=(12,12))
# for nCase,Case in enumerate(RunList):
#     plt.subplot(3,3,nCase+1)
#     plt.plot(sounding[nCase,:],interpZ)
#     plt.axhline(y=zi[nCase], ls='--',label='zi=%d' %zi[nCase])
#     plt.title('%s: PHI=%.2f' %(Case,Phi[nCase]))
#     plt.gca().set(xlabel='Theta [K]',ylabel='z [m]',ylim =[0,2500])
#     plt.legend()
# plt.tight_layout()
# # plt.show()
# plt.savefig(plume.figdir + 'CWIzCL/ziTEST.pdf')
