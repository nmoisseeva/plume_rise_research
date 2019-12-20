# nmoisseeva@eoas.ubc.ca
# December 2019

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

g = 9.81

r = np.empty((runCnt)) * np.nan
Ua = np.empty((runCnt)) * np.nan
H = np.empty((runCnt)) * np.nan
zi = np.empty((runCnt)) * np.nan
Ti = np.empty((runCnt)) * np.nan
Gamma = np.empty((runCnt)) * np.nan
zCL = np.empty((runCnt)) * np.nan
wStar = np.empty((runCnt)) * np.nan
zCLend = np.empty((runCnt)) * np.nan

for nCase,Case in enumerate(RunList):
    if Case in plume.exclude_runs:
        continue
    print('Examining case: %s ' %Case)

    #----------get preignition temperature profile-----------------
    cspath = plume.wrfdir + 'interp/wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    #save initial temperature prfile
    profpathT = plume.wrfdir + 'interp/profT0' + Case + '.npy'
    profileT = np.mean(csdict['temp'][0,:,:],1)
    np.save(profpathT,profileT)

    profpathU = plume.wrfdir + 'interp/profU0' + Case + '.npy'
    profileU = np.mean(csdict['u'][0,:,:],1)
    np.save(profpathU,profileU)

    #----------check for interpolated data----------------------------
    avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'

    if os.path.isfile(avepath):
        print('Averaged data found at: %s' %avepath)
        avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
    else:
        sys.exit('ERROR: no averaged data found - run prep_plumes.py via submit_interp.sh first on Cedar!')

    pm = ma.masked_where(avedict['pm25'] <= 30, avedict['pm25'] )
    w = ma.masked_where(avedict['pm25'] <= 30, avedict['w'] )
    # pm = ma.masked_where(avedict['pm25'] <= avedict['pm25'].max()*0.001, avedict['pm25'] )
    # w = ma.masked_where(avedict['pm25'] <= avedict['pm25'].max()*0.001, avedict['w'] )

    PMmaxVidx = pm.argmax(0)
    xmax,ymax = np.nanargmax(PMmaxVidx), np.nanmax(PMmaxVidx)
    # centerline = ma.masked_where(plume.lvl[PMmaxVidx] == 0, plume.lvl[PMmaxVidx])
    centerline = plume.lvl[PMmaxVidx]
    zCLend[nCase] = np.mean(centerline[centerline > 0][-5:])

    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')


    ignited = np.array([i for i in avedict['ghfx'] if i > 2])
    r[nCase] = len(ignited) * plume.dx
    H[nCase] = np.mean(ignited) * 1000 / ( 1.2 * 1005)         #get heat flux


    #calculate injection height based on temperature anomaly
    diffT = (avedict['temp'].T-T0).T
    Tctr = np.array([diffT[i,ni] for ni, i in enumerate(PMmaxVidx)])
    Wctr = np.array([w[i,ni] for ni, i in enumerate(PMmaxVidx)])
    intXnan = np.argwhere(np.diff(np.sign(Tctr))).flatten()             #find where the interesepts are
    intX = intXnan[np.isfinite(Tctr[intXnan])]                          #get only the values in the plume
    intX = intX[intX > plume.wi]                                        #get rid of everything ahead of fireline
    intX = intX[PMmaxVidx[intX] > 0]                                    #remove any near-surface kinks
    if len(intX) < 2:
        sliceX = np.nanargmin(Wctr[plume.wi:]) + plume.wi
    else:
        sliceX = intX[1]
    sliceZ = PMmaxVidx[sliceX]

    zCL[nCase] = plume.lvl[sliceZ]
    if zCL[nCase] < 500 :
        print('\033[93m' + 'Plume zCL: %d m' %zCL[nCase] + '\033[0m')
        print(zCLend[nCase])
    else:
        print('Plume zCL: %d m' %zCL[nCase])
        print(zCLend[nCase])


    #BL characteristics -------------------------------------------------
    dT = T0[1:]-T0[0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 3                                                      #surface layer
    Ti[nCase] = T0[si+1]                                        #characteristic BL temperature
    zi_idx = np.argmax(gradT[si:]) + si                         #vertical level index of BL top
    Gamma[nCase] = np.mean(dT[zi_idx:zi_idx+20])/plume.dz       #inversion strength
    print('inversion std: %0.2f' %np.std(dT[zi_idx:zi_idx+20]))
    zi[nCase] = plume.dz * zi_idx                               #BL height
    Ua[nCase] = np.mean(U0[si:10])                               #ambient BL wind

    #analysis with w*
    wStar[nCase] = (g*zi[nCase]*H[nCase]*Gamma[nCase]/(Ti[nCase]))**(1/4.)





wStar_bar = wStar / Ua
zCL_bar = zCLend / r

fig = plt.figure(figsize=(12,4))
plt.suptitle('DIMENSIONLESS ANALYSIS')
plt.subplot(1,3,1)
plt.scatter(wStar_bar, zCL_bar, c=plume.read_tag('F',RunList),cmap=plt.cm.tab20)
# plt.xlabel(r'$\frac{H}{\Gamma * z_i * w_*}$')
# plt.ylabel(r'$\frac{z_{CL}}{z_i}$ ')
plt.colorbar()


plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'DimAnalysis_Wf.pdf' )
plt.close()
