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
from scipy.stats import linregress

#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time
zstep = 20              #height interpolation step

#=================end of input===============

RunList =   [i for i in plume.tag if i not in plume.exclude_bad]
# RunList =   ['W5F4R5TE','W3F7R7T','W4F4R1']
runCnt = len(RunList)                           #count number of cases

dropoff = np.empty((runCnt)) * np.nan                #BL height (m)
Phi = np.empty((runCnt)) * np.nan
wf = np.empty((runCnt)) * np.nan
span = np.empty((runCnt)) * np.nan
#======================repeat main analysis for all runs first===================
#loop through all LES cases
plumes = np.load('plumeData.npy',allow_pickle=True).item()
for nCase,Case in enumerate(RunList):
    print('...case: %s' %Case)
    #create an interpolated profile of temperature
    if Case[-1:]=='E':
        continue
    elif Case[-1:]=='T':
        levels = plume.lvltall
        pmlvl=levels
    else:
        levels=plume.lvl
        pmlvl=levels

    runNumber = np.where(np.array(plumes['RunList'])==Case)[0][0]
    zCL =  plumes['zCL'][runNumber]   #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(pmlvl - zCL))
    zi =  plumes['zi'][runNumber]
    ziidx = np.argmin(abs(pmlvl - zi))

    if zCL < zi + zstep:
        continue
    else:
        #load 3D fields
        endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
        enddict = np.load(endpath, allow_pickle=True).item()
        dimZ,dimY,dimX = np.shape(enddict['PM25'][:,:,:])
        lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
        crosszCL = lateralPM[zCLidx,:]
        normMax = np.nanmax(crosszCL)
        NlateralPM = lateralPM/normMax
        NcrosszCL = crosszCL/normMax
        yaxis = np.arange(0, int(dimY * plume.dy), int(plume.dy))



        fig = plt.figure(figsize=(6,5))
        plt.subplot(211)
        plt.title('(a) ALONG-WIND TOTAL SMOKE')
        im = plt.imshow(NlateralPM, origin='lower', extent=[0,yaxis[-1],0,pmlvl[-1]],cmap=plt.cm.cubehelix_r,vmin=0,vmax = 1)
        plt.colorbar(im, label=r'normalized concentration')
        plt.axhline(y = zCL,ls='--', c='dimgrey',label=r'$z_{CL}$' )
        plt.gca().set(ylim=[0,pmlvl[-1]],aspect='equal',ylabel='height [m]',xlabel='distance [m]')
        plt.legend()
        plt.subplot(212)
        plt.title(r'(b) CROSSSECTION AT $z_{CL}$')
        plt.plot(yaxis,NcrosszCL)
        plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
        plt.tight_layout()
        # plt.show()
        plt.savefig(plume.figdir + 'lateral/lat%s.pdf' %Case )
        plt.close()

        #get plume edges
        edge1 = np.nanargmin(abs(NcrosszCL[:int(dimY/2)] - 0.21)) * plume.dy
        edge2 = (np.nanargmin(abs(NcrosszCL[int(dimY/2):] - 0.21)) + int(dimY/2))* plume.dy
        dist1 = 4000 - edge1
        dist2 = edge2 - 6000
        dropoff[nCase] = np.mean([dist1,dist2])
        Phi[nCase] = plumes['Phi'][runNumber]
        # wf[nCase] = Phi[nCase]/(plumes['zi'][runNumber] * plumes['Omega'][runNumber])
        wf[nCase] = (9.81 * Phi[nCase] * (plumes['zCL'][runNumber] - 0.75 * plumes['zi'][runNumber] )/ (plumes['zi'][runNumber] * plumes['thetaS'][runNumber]))**(1/3.)

        span[nCase] = np.mean([edge2, edge1])
plt.figure()
xaxis = np.arange(0,np.nanmax(Phi))
plt.title('PLUME WIDENING')
plt.scatter(Phi,dropoff,c=plume.read_tag('F',RunList))
plt.plot(xaxis, 1000+xaxis*(2500)/(140000),c='C1',ls='--',label='linear fit')
plt.gca().set(xlabel=r'fireline intensity [K m$^2$s$^{-1}$]',ylabel='plume widening [m]')
plt.legend()
plt.savefig(plume.figdir + 'lateral/Widening.pdf')
plt.show()

plt.figure()
xaxis = np.arange(0,np.nanmax(Phi))
plt.title('PLUME WIDENING')
plt.scatter(wf,span,c=plume.read_tag('F',RunList))
plt.gca().set(xlabel=r'$w_f$',ylabel='plume span [m]')
plt.legend()
fitWf = linregress(wf[np.isfinite(wf)], span[np.isfinite(span)])

#======================compare conditions:Legth===================
clr = ['C0','C1']
#length tests
RunList =   ['W4F7R4L1','W4F7R4L4']
RunLbl =   ['1 km','4 km']
edgeDist = [113,75]
runCnt = len(RunList)                           #count number of cases

#loop through all LES cases
plt.figure(figsize = (6,10))
plt.subplot(3,1,1)
plt.title(r'(a) FIRELINE LENGTH EFFECT')
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
    PhiLoc = plumes['Phi'][runNumber]


    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))

    #modelling class
    sigma = 1000+PhiLoc*(2500)/(140000)
    sigmaG = (500*wf[runNumber] + 2000)
    modelledCS = np.empty_like(NcrosszCL) * np.nan
    s, f = edgeDist[nCase], (dimY-edgeDist[nCase])
    modelledCS[s:f] =0.9
    modelledCS[:s] =0.9* np.exp(-0.5*((yaxis[:s] - s*plume.dy)/sigma)**2)
    modelledCS[f:] =0.9* np.exp(-0.5*((yaxis[f:] - f*plume.dy)/sigma)**2)
    GaussianCS = np.exp((-0.5*((yaxis - 5000)/sigmaG)**2))
    plt.plot(yaxis,GaussianCS, ls='--',c=clr[nCase])

    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])
    # plt.plot(yaxis,modelledCS, ls=':',c=clr[nCase])
    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()


#======================compare conditions:Wind===================

#length tests
RunList =   ['W3F7R2','W12F7R2']
RunLbl =   ['3 m/s','12 m/s']
runCnt = len(RunList)                           #count number of cases
edgeDist = 100

#loop through all LES cases
plt.subplot(3,1,2)
plt.title(r'(b) AMBIENT WIND EFFECT')
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
    PhiLoc = plumes['Phi'][runNumber]

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))

    #modelling class
    sigma = 1000+PhiLoc*(2500)/(140000)
    sigmaG = (500*wf[runNumber] + 2000)
    modelledCS = np.empty_like(NcrosszCL) * np.nan
    s, f = edgeDist, (dimY-edgeDist)
    modelledCS[s:f] = 0.9
    modelledCS[:s] = 0.9*np.exp(-0.5*((yaxis[:s] - s*plume.dy)/sigma)**2)
    modelledCS[f:] = 0.9*np.exp(-0.5*((yaxis[f:] - f*plume.dy)/sigma)**2)
    GaussianCS = np.exp((-0.5*((yaxis - 5000)/sigmaG)**2))
    plt.plot(yaxis,GaussianCS, ls='--',c=clr[nCase])

    # plt.plot(yaxis,modelledCS, ls=':',c=clr[nCase])
    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])


    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()


#======================compare conditions:Intensity===================

#length tests
RunList =   ['W5F1R1','W5F6R1']
RunLbl =   ['low intensity','high intensity']
runCnt = len(RunList)                           #count number of cases

#loop through all LES cases
plt.subplot(3,1,3)
plt.title(r'(c) FIRELINE HEAT EFFECT')
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
    PhiLoc = plumes['Phi'][runNumber]

    #load 3D fields
    endpath = plume.wrfdir+'interp/end/wrfend_' + Case + '.npy'
    enddict = np.load(endpath, allow_pickle=True).item()
    lateralPM = np.nansum(enddict['PM25'][:,:,:],2)
    crosszCL = lateralPM[zCLidx,:]
    normMax = np.nanmax(crosszCL)
    NlateralPM = lateralPM/normMax
    NcrosszCL = crosszCL/normMax
    yaxis = np.arange(0, int(np.shape(lateralPM)[1] * plume.dy), int(plume.dy))

    #modelling class
    sigma = 1000+PhiLoc*(2500)/(140000)
    sigmaG = (500*wf[runNumber] + 2000)
    modelledCS = np.empty_like(NcrosszCL) * np.nan
    s, f = edgeDist, (dimY-edgeDist)
    modelledCS[s:f] = 0.9
    modelledCS[:s] = 0.9*np.exp(-0.5*((yaxis[:s] - s*plume.dy)/sigma)**2)
    modelledCS[f:] = 0.9*np.exp(-0.5*((yaxis[f:] - f*plume.dy)/sigma)**2)
    GaussianCS = np.exp((-0.5*((yaxis - 5000)/sigmaG)**2))
    plt.plot(yaxis,GaussianCS, ls='--',c=clr[nCase])

    plt.plot(yaxis,NcrosszCL,label=RunLbl[nCase])
    # plt.plot(yaxis,modelledCS, ls=':',c=clr[nCase])
    plt.gca().set(ylabel=r'normalized concentration',xlabel='distance [m]')
    plt.legend()

# plt.show()
plt.tight_layout()
plt.savefig(plume.figdir + 'lateral/CompareEffects.pdf')
plt.close()
