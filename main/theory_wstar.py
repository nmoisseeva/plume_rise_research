# nmoisseeva@eoas.ubc.ca
# October 2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time


#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
runCnt = len(RunList)


radius = np.empty((runCnt,len(plume.lvl))) * np.nan
wCL = np.empty((runCnt,len(plume.lvl))) * np.nan
B = np.empty((runCnt,len(plume.lvl))) * np.nan
Uz = np.empty((runCnt,len(plume.lvl))) * np.nan

for nCase,Case in enumerate(RunList):
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)

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

    #mask plume as being at last 30ppm---------------------------------
    # pm = ma.masked_where(avedict['pm25'] <= 30, avedict['pm25'] )
    # w = ma.masked_where(avedict['pm25'] <= 30, avedict['w'] )
    pm = avedict['pm25']
    w = avedict['w']

    PMmaxVidx = pm.argmax(0)
    xmax,ymax = np.nanargmax(PMmaxVidx), np.nanmax(PMmaxVidx)
    centerline = ma.masked_where(plume.lvl[PMmaxVidx] == 0, plume.lvl[PMmaxVidx])

    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')
    Uz[nCase,:] = U0

    #charatertistics of plume temperature anomalyies---------------------
    # diffT = ma.masked_where(avedict['pm25'] <= 30, (avedict['temp'].T-T0).T)           #get temperature anomaly
    diffT = (avedict['temp'].T-T0).T
    Tctr = np.array([diffT[i,ni] for ni, i in enumerate(PMmaxVidx)])
    Wctr = np.array([w[i,ni] for ni, i in enumerate(PMmaxVidx)])
    Uctr = np.array([avedict['u'][i,ni] for ni, i in enumerate(PMmaxVidx)])

    # #calculate total flux from ground and fire
    ignited = np.array([i for i in avedict['ghfx'] if i > 2])
    r = len(ignited) * plume.dx

    #charageteristic vertical velocity--------------------------------------------
    g = 9.81
    Fflux = np.mean(ignited) * 1000 / ( 1.2 * 1005)         #get heat flux
    # Wf = (g*BLdict['zi'][nCase]*Fflux/BLdict['Ti'][nCase])**(1./3)                     #characteristic velocity based on Deodorff

    z1=np.sqrt( (Fflux * g * r**3) / (T0[0] * U0[3]**3) )

    # A[nCase] = g* parcelHeat[nCase]/ (zi*Ti[nCase])           #some sort of non-dimensional variable


    #get profile location---------------------------------------------------------
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


    #test entrainment--------------------------------------------------------------
    #test whether plume width entrainment is a function of vertcal VELOCITY

    #get plume width at each level
    for nLvl in range(1,len(plume.lvl)):
        hslice = pm[nLvl,:]
        xCL = np.nanargmax(hslice)
        edge = hslice.max() * 0.1              #cutoff at 10% of maximum at that level
        edgeLidx = np.nanargmin(abs(hslice - edge)[0:xCL])
        edgeRidx = np.nanargmin(abs(hslice - edge)[xCL:]) + xCL
        if nLvl == sliceZ:
            break
        radius[nCase,nLvl] = (edgeRidx - edgeLidx) * plume.dx
        wCL[nCase,nLvl] = np.nanmax(w[nLvl,:])
        tempCL = np.nanargmax(diffT[nLvl,:])
        B[nCase,nLvl] = g* tempCL / T0[nLvl]

    fig = plt.figure(figsize=(10,10))
    plt.suptitle('PLUME RADIUS: %s' %Case)
    plt.subplot(2,2,1)
    plt.scatter(wCL[nCase,:nLvl],radius[nCase,:nLvl])
    plt.ylabel('radius [m]')
    plt.xlabel('w [m/s]')

    plt.subplot(2,2,2)
    plt.scatter(1/wCL[nCase,:nLvl]**2,radius[nCase,:nLvl])
    plt.ylabel('radius [m]')
    plt.xlabel( '$1/w^2$  $[s^2/m^2]$')

    plt.subplot(2,2,3)
    plt.scatter(B[nCase,:nLvl]/wCL[nCase,:nLvl],radius[nCase,:nLvl])
    plt.ylabel('radius [m]')
    plt.xlabel('$B/w^2$ $[m^{-1}]$')

    plt.subplot(2,2,4)
    plt.scatter(U0[:nLvl],radius[nCase,:nLvl])
    plt.ylabel('radius [m]')
    plt.xlabel('Uz [m/s]')

    plt.subplots_adjust(top=0.85)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(plume.figdir + 'entrainment/e%s.pdf' %Case)
    plt.close()

fig = plt.figure(figsize=(10,10))
plt.suptitle('PLUME RADIUS: %s' %Case)
plt.subplot(2,2,1)
plt.scatter(wCL[:,:nLvl].ravel(),radius[:,:nLvl].ravel())
plt.ylabel('radius [m]')
plt.xlabel('w [m/s]')

plt.subplot(2,2,2)
plt.scatter((1/wCL[:,:nLvl]**2).ravel(),radius[:,:nLvl].ravel())
plt.ylabel('radius [m]')
plt.xlabel( '$1/w^2$  $[s^2/m^2]$')

plt.subplot(2,2,3)
plt.scatter(B[:,:nLvl].ravel()/wCL[:,:nLvl].ravel(),radius[:,:nLvl].ravel())
plt.ylabel('radius [m]')
plt.xlabel('$B/w^2$ $[m^{-1}]$')

plt.subplot(2,2,4)
plt.scatter(Uz[:,:nLvl].ravel()/wCL[:,:nLvl].ravel(),radius[:,:nLvl].ravel())
plt.ylabel('radius [m]')
plt.xlabel('Uz [m/s]$')


# plt.subplots_adjust(top=0.85)
# plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.show()
# #plt.savefig(plume.figdir + 'entrainment/e%s.pdf' %Case)
