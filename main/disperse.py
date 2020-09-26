# September 2020
#nmoisseeva@eoas.ubc.ca
#This code laterally disperses downwind-averaged concentrations

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from matplotlib import gridspec
from scipy.interpolate import interp1d


#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time
#=================end of input===============

# RunList =   [i for i in plume.tag if i not in plume.exclude_bad]
RunList =   ['W5F4R5TE','W3F7R7T','W4F4R1']
runCnt = len(RunList)                           #count number of cases


#======================repeat main analysis for all runs first===================
#loop through all LES cases
plumes = np.load('plumeData.npy',allow_pickle=True).item()
for nCase,Case in enumerate(RunList):
    print('...case: %s' %Case)
    #create an interpolated profile of temperature
    if Case[-1:]=='T' or Case[-1:]=='E':
        levels = plume.lvltall
        if Case[-1:]=='E':
            pmlvl = np.arange(0,4001,40)        #extra tall domain for the concentrations only
        else:
            pmlvl=levels
    else:
        levels=plume.lvl
        pmlvl=levels

    runNumber = np.where(np.array(plumes['RunList'])==Case)[0][0]
    zCL =  plumes['zCL'][runNumber]   #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(pmlvl - zCL))
    zi =  plumes['zi'][runNumber]
    ziidx = np.argmin(abs(pmlvl - zi))

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))


    
    # fig = plt.figure(figsize=(6,5))
    # plt.subplot(211)
    # plt.title('LATERAL MEAN')
    # im = plt.imshow(NlateralPM, origin='lower', extent=[0,yaxis[-1],0,pmlvl[-1]],cmap=plt.cm.cubehelix_r,vmin=0,vmax = 1)
    # plt.colorbar(im, label=r'normalized concentration')
    # plt.axhline(y = zCL,ls='--', c='dimgrey',label=r'$z_{CL}$' )
    # plt.gca().set(ylim=[0,pmlvl[-1]],aspect='equal',ylabel='height [m]',xlabel='distance [m]')
    # plt.legend()
    # plt.subplot(212)
    # plt.title(r'CROSSSECTION AT $z_{CL}$')
    # plt.plot(yaxis,NcrosszCL)
    # plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(plume.figdir + 'lateral/lat%s.pdf' %Case )
    # plt.close()

#======================compare conditions:Legth===================

#length tests
RunList =   ['W4F7R4L1','W4F7R4L4']
RunLbl =   ['1 km','4 km']
runCnt = len(RunList)                           #count number of cases

#loop through all LES cases
plt.figure(figsize = (6,10))
plt.subplot(3,1,1)
plt.title(r'LENGTH EFFECT')
for nCase,Case in enumerate(RunList):
    print('...case: %s' %Case)
    #create an interpolated profile of temperature
    if Case[-1:]=='T' or Case[-1:]=='E':
        levels = plume.lvltall
        if Case[-1:]=='E':
            pmlvl = np.arange(0,4001,40)        #extra tall domain for the concentrations only
        else:
            pmlvl=levels
    else:
        levels=plume.lvl
        pmlvl=levels

    runNumber = np.where(np.array(plumes['RunList'])==Case)[0][0]
    zCL =  plumes['zCL'][runNumber]   #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(pmlvl - zCL))
    zi =  plumes['zi'][runNumber]
    ziidx = np.argmin(abs(pmlvl - zi))

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))
    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])
    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()


#======================compare conditions:Wind===================

#length tests
RunList =   ['W3F7R2','W12F7R2']
RunLbl =   ['3 m/s','12 m/s km']
runCnt = len(RunList)                           #count number of cases

#loop through all LES cases
plt.subplot(3,1,2)
plt.title(r'AMBIENT WIND EFFECT')
for nCase,Case in enumerate(RunList):
    print('...case: %s' %Case)
    #create an interpolated profile of temperature
    if Case[-1:]=='T' or Case[-1:]=='E':
        levels = plume.lvltall
        if Case[-1:]=='E':
            pmlvl = np.arange(0,4001,40)        #extra tall domain for the concentrations only
        else:
            pmlvl=levels
    else:
        levels=plume.lvl
        pmlvl=levels

    runNumber = np.where(np.array(plumes['RunList'])==Case)[0][0]
    zCL =  plumes['zCL'][runNumber]   #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(pmlvl - zCL))
    zi =  plumes['zi'][runNumber]
    ziidx = np.argmin(abs(pmlvl - zi))

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))
    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])
    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()



#======================compare conditions:Intensity===================

#length tests
RunList =   ['W5F1R2','W5F6R2']
RunLbl =   ['low intensity','high intensity']
runCnt = len(RunList)                           #count number of cases

#loop through all LES cases
plt.subplot(3,1,3)
plt.title(r'FIRELINE HEAT EFFECT')
for nCase,Case in enumerate(RunList):
    print('...case: %s' %Case)
    #create an interpolated profile of temperature
    if Case[-1:]=='T' or Case[-1:]=='E':
        levels = plume.lvltall
        if Case[-1:]=='E':
            pmlvl = np.arange(0,4001,40)        #extra tall domain for the concentrations only
        else:
            pmlvl=levels
    else:
        levels=plume.lvl
        pmlvl=levels

    runNumber = np.where(np.array(plumes['RunList'])==Case)[0][0]
    zCL =  plumes['zCL'][runNumber]   #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(pmlvl - zCL))
    zi =  plumes['zi'][runNumber]
    ziidx = np.argmin(abs(pmlvl - zi))

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))
    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])
    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()

# plt.show()
plt.tight_layout()
plt.savefig(plume.figdir + 'lateral/CompareEffects.pdf')
plt.close()
