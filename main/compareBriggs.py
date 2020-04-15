# April 2020
#nmoisseeva@eoas.ubc.ca
#This code applies FEPS version of Briggs to LES data: https://www.fs.fed.us/pnw/fera/feps/FEPS_users_guide.pdf
# It also partitions LES runs into model and test sets and applies the injection height parameterization
# Plotting shows model sensitivity and error distributions

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.optimize import fsolve
from matplotlib import gridspec

#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

testPortion = 0.2       #portion of data to reserve for testing the model
trials = 10             #number of times to rerun the model

g = 9.81                #gravity constant
si = 3                  #number of levels to skip from the bottom (exlude surface layer)
#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
# RunList = ['W5F3R3']
runCnt = len(RunList)


#saved variables
r = np.empty((runCnt)) * np.nan                 #fireline depth
Ua = np.empty((runCnt)) * np.nan                #ambient wind
zi = np.empty((runCnt)) * np.nan                #BL height
zCL = np.empty((runCnt)) * np.nan               #centerline height
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature
Phi = np.empty((runCnt)) * np.nan               #cumulative fire heat
BL = np.empty((runCnt,len(plume.lvl)-1))        #storage for BL lapse rates
FlaggedCases = []                               #for storage of indecies of anomalous runs

#Briggs variables
Wi = np.empty((runCnt)) * np.nan                #transport wind speed
s = np.empty((runCnt)) * np.nan                 #stability parameter
Fi = np.empty((runCnt)) * np.nan                #Briggs buoyancy factor
BriggsH = np.empty((runCnt)) * np.nan           #Briggs plume rise

#======================repeat main LES analysis for all runs first===================

for nCase,Case in enumerate(RunList):
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)

    csdict = plume.load_CS_prep_Profiles(Case)
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    #mask plume with cutoff value---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])
    zi[nCase] = plume.get_zi(T0)
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= plume.PMcutoff, csdict['pm25'][-1,:,:] )

    #locate centerline
    ctrZidx = pm.argmax(0)
    ctrXidx = pm.argmax(1)
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])
    for nZ in range(dimZ):
        if pmCtr[ctrXidx[nZ]] < pm[nZ,ctrXidx[nZ]]:
            pmCtr[ctrXidx[nZ]] = pm[nZ,ctrXidx[nZ]]
    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
    centerline = ma.masked_where(plume.lvl[ctrZidx] == 0, plume.lvl[ctrZidx])
    smoothCenterline = savgol_filter(centerline, 31, 3)         # window size 31, polynomial order 3

    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 101, 3)                     # window size 101, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.05 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]
    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.mean(stablePM,1)
    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)

    #define source (r and H)------------------------
    #raduis using full 2D average -
    masked_flux = ma.masked_less_equal(csdict['ghfx2D'], 0.5)
    cs_flux = np.nanmean(masked_flux,1)                         #get cross section for each timestep
    fire = []
    xmax = np.argmax(cs_flux,axis=1)                            #get max for each timestep
    for nP, pt in enumerate(xmax[plume.ign_over:]):             #excludes steps containing ignition
        subset = cs_flux[plume.ign_over+nP,pt-plume.wi:pt+plume.wf]         #set averaging window around a maximum
        fire.append(subset)
    burning = np.sum(np.sum(csdict['ghfx2D'][plume.ign_over:,:,:],1),1)     #get total heat flux from entire fire area

    meanFire = np.nanmean(fire,0)                               #geta verage cross-section
    ignited = np.array([i for i in meanFire if i > 0.5])        #extract ignited grids (ignited: >0.5 kW/m2)
    r[nCase] = len(ignited) * plume.dx                          #get cross-length of fire
    Phi[nCase] = np.trapz(ignited, dx = plume.dx) * 1000 / ( 1.2 * 1005)    #calculate intgrated fireline using kinetic flux

    #injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is region of stable centerline and concentration
    zCLidx = int(np.mean(ctrZidx[1:][stablePMmask]))            #index of model level closest to centerline height
    dT = T0[1:]-T0[0:-1]                                        #get "lapse rate" profileT (units of K ONLY!)
    Omega[nCase] = np.trapz(dT[si+1:zCLidx], dx = plume.dz)     #calculate the integrated profile term
    BL[nCase,:] = dT                                            #store "lapse rate" profile
    Ua[nCase] = np.mean(U0[si:zCLidx])                          #store mean wind


    #----------------------FEPS BRIGGS PORTION-------------------------
    #get transport wind speed (Eq 34)
    Wi[nCase] = max([5, Ua[nCase]])                             #assuming no previous hour speed/same as current

    #get Briggs stability parameter (Eq 35)
    Ti = T0[si+1]                                               #get temperature
    ziidx = np.argmin(abs(plume.lvl - zi[nCase]))               #get index of BL
    dTBriggs = (T0[-1] - T0[ziidx]) / (plume.lvl[-1]- plume.lvl[ziidx])     #calculate lapse rate
    s[nCase] = g * dTBriggs / Ti                                #save stability parameter

    #get heat release (Eq 39)
    fire2D = masked_flux[40,:,:]                                #get the fire 2D flux half way through the burn
    fireKW = np.nansum(fire2D)*(plume.dx**2)                      #get total fire KW
    Qi = fireKW * 3412.142                                      #convert to BTU/h

    #get Briggs Buoyancy factor (Eq 39)
    Fi[nCase] = Qi * 0.00000258                                 #units of M4/sec3

    #do Briggs Calculation (Eq 43 and 44)
    condition1 = Fi[nCase]**(0.6) * 38.7 / Wi[nCase]
    condition2 = 2.4 * ((Fi[nCase] / (Wi[nCase] * s[nCase]))**0.333) * 38.7 / Wi[nCase]
    if s[nCase] >= 0:
        BriggsH[nCase] = min([condition1,condition2])
    else:
        BriggsH[nCase] = condition1
    print(BriggsH[nCase])

#======================plot Briggs error===================

plt.figure()
plt.title('LES EVALUATION: FEPS BRIGGS')
ax = plt.gca()
plt.scatter(zCL, BriggsH)
ax.set(xlabel='LES $z_{CL}$ [m]', ylabel='FEPS Briggs $\Delta H_{max}$ [m]')
plt.savefig(plume.figdir + 'injectionModel/BriggsPerformance.pdf')
plt.show()
