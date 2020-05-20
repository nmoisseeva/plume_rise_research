#April 2020
#nmoisseeva@eoas.ubc.ca

#This code explores relationship bewtween fire frequency Nf and plume tilt
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
Phi = np.empty((runCnt)) * np.nan            #cumulative fire intensity
tilt = np.empty((runCnt)) * np.nan          #plume tilt
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
    Phi[nCase] = np.trapz(ignited, dx = plume.dx)

    #mask plume with cutoff value---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])
    ave_pm = np.nanmean(csdict['pm25'][-20:,:,:],0)
    pm = ma.masked_where(ave_pm <= plume.PMcutoff, ave_pm )
    zi[nCase] = plume.get_zi(T0)

    #locate centerline
    ctrZidx = pm.argmax(0)
    ctrZidx[:fmax[-1]] = 0
    ctrXidx = pm.argmax(1)

    #define downwind distribution
    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
    tilt[nCase] = ymax/xmax

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

    # #get final profile weights by calculating proportional area for each level
    # subArea = stableProfile * plume.dx
    # totalArea = sum(subArea.data)
    # weights = subArea/totalArea


    #dimensional analysis variables ---------------------------
    dT = T0[1:]-T0[0:-1]
    Ua[nCase] = np.mean(U0[si:zCLidx])
    Omega[nCase] = np.trapz(dT[si+1:zCLidx], dx = plume.dz)

##now plot zCl as a function of w*
Nf = np.sqrt(g*Phi/(Omega *zi* Ua))



#plot "fire frequency" vs tilt
slope, intercept, r_value, p_value, std_err = linregress(Nf,tilt)
print('Sum of residuals: %0.2f' %r_value)

fig = plt.figure(figsize=(12,6))
# plt.suptitle('zCL=FCN(W*): R = %.2f' %r_value)
plt.subplot(1,2,1)
ax1=plt.gca()
plt.scatter(Nf, tilt,  c=Ua, cmap=plt.cm.jet)
# plt.plot(Nf, intercept + slope*Nf, c='grey', label='R = %.2f' %r_value)
plt.colorbar(label='Ua wind speed [m/s]')
ax1.set(xlabel='Nf',ylabel='plume tilt')
# plt.legend(loc='upper left')
plt.subplot(1,2,2)
ax2=plt.gca()
plt.scatter(Nf, tilt,  c=zi, cmap=plt.cm.copper)
plt.colorbar(label='zi [m]')
ax2.set(xlabel='Nf',ylabel='plume tilt')
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig(plume.figdir + 'zCl_wStar.pdf' )
plt.show()
plt.close()
