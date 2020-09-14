# March 2020
#nmoisseeva@eoas.ubc.ca
# This code partitions LES runs into model and test sets and applies the injection height parameterization
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
from scipy.interpolate import interp1d
from scipy.stats import norm


#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time

g = 9.81                #gravity constant
zstep = 20              #height interpolation step
BLfrac = 0.75           #fraction of BL height to set zs at
#=================end of input===============

RunList =   [i for i in plume.tag if i not in plume.exclude_bad]
runCnt = len(RunList)                           #count number of cases

#set up interpolated vertical profile with 5m vertical step
interpZ = np.arange(0, 4001, zstep)

#storage for variables
zi = np.empty((runCnt)) * np.nan                #BL height (m)
zCL = np.empty((runCnt)) * np.nan               #smoke injection height (m)
Phi = np.empty((runCnt)) * np.nan               #cumulative fire heat  (K m^2/s)
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature (Km) - the questionable denominator term
thetaS = np.empty((runCnt)) * np.nan
thetaCL = np.empty((runCnt)) * np.nan
profile = np.empty((runCnt,len(interpZ))) * np.nan
quartiles = np.empty((runCnt,len(interpZ),2)) * np.nan


sounding = np.empty((runCnt,len(interpZ))) * np.nan         #storage for interpolated soundings
gradT0interp = np.empty((runCnt,len(interpZ)-1)) * np.nan   #storage for temperature gradient

#======================repeat main analysis for all runs first===================
#loop through all LES cases
for nCase,Case in enumerate(RunList):
    csdict = plume.load_CS_prep_Profiles(Case)      #load data cross-wind integrated smoke and all other data
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')    #load initial temperature profile
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')    #load intial wind profile

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
    centerline = ma.masked_where(pmlvl[ctrZidx] == 0, pmlvl[ctrZidx])               #make sure centerline is only calculated inside the plume
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
                            centerline[nX] < pmlvl[-1]-200 and \
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

    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.median(stablePM,1)

    profile[nCase,:] = interp1d(pmlvl,stableProfile,fill_value='extrapolate')(interpZ)

    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)
    quartiles[nCase,:,0] = interp1d(pmlvl,pmQ1,fill_value='extrapolate')(interpZ)
    quartiles[nCase,:,1] = interp1d(pmlvl,pmQ3,fill_value='extrapolate')(interpZ)

    #define heat source * to include larger fires use padding ------------------------
    masked_flux_padded = ma.masked_less_equal(np.pad(csdict['ghfx2D'],((0,0),(0,0),(100,0)), 'constant',constant_values=0),1)
    cs_flux_padded = np.nanmean(masked_flux_padded,1)
    fire_padded = []
    fxmax_padded = np.argmax(cs_flux_padded,axis=1)
    for nP, pt in enumerate(fxmax_padded[plume.ign_over:]):            #excludes steps containing ignition
        subset_padded = cs_flux_padded[plume.ign_over+nP,pt-100:pt+plume.wf]     #set averaging window around a maximum
        fire_padded.append(subset_padded)

    meanFire_padded = np.nanmean(fire_padded,0)                               #calculate mean fire cross section
    ignited_padded = np.array([i for i in meanFire_padded if i > 0.5])        #consider only cells that have heat flux about 500 W/m2
    Phi[nCase] = np.trapz(ignited_padded, dx = plume.dx) * 1000 / ( 1.2 * 1005)    #calculate Phi by integrating kinematic heat flux along x (Km2/s)

    #calculate injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
    ziidx = np.argmin(abs(interpZ - zi[nCase]))
    dT = (T0[1:]-T0[0:-1])/plume.dz                                        #calculate potential temperature change (K)

    zsidx = np.argmin(abs(interpZ - zi[nCase]*BLfrac))
    sounding[nCase,:] = T0interp
    dTinterp = (T0interp[1:] - T0interp[0:-1])/zstep
    gradT0interp[nCase,:] = dTinterp
    Omega[nCase] = np.trapz(dTinterp[zsidx:zCLidx],dx=zstep)
    thetaS[nCase] = sounding[nCase,zsidx]
    thetaCL[nCase] = sounding[nCase,zCLidx]

#================save data to avoid rerunning=======
plume_dict = {'RunList':RunList, 'zCL':zCL, 'zi':zi, 'profile':profile, 'sounding':sounding, \
                        'Phi':Phi, 'thetaS':thetaS, 'thetaCL':thetaCL, 'Omega':Omega}
np.save('plumeData.npy', plume_dict)

#===========iterative solution===============
zS = zi*BLfrac
mf, bf = 0.9243, 115.059
C = 1.0087

zCLerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
zCLmodel = np.empty((runCnt)) * np.nan
for nCase,Case in enumerate(RunList):
    BLidx = np.nanargmin(abs(interpZ - BLfrac*zi[nCase]))

    toSolveCase = lambda z : z  - bf - mf*(zS[nCase] + \
                    C/(np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
                    (g*Phi[nCase]*(z-zi[nCase]*BLfrac)/(thetaS[nCase] * zi[nCase]))**(1/3.))

    zCLmodel[nCase] = fsolve(toSolveCase, zi[nCase],factor=0.1)             #solve
    zCLerror[nCase] =  zCL[nCase]  - zCLmodel[nCase]

plt.scatter(zCLmodel,zCL)
plt.gca().set(aspect='equal')
#=============find penetrative plumes========
ABLidx = np.where(zi+zstep > zCL)[0]
PENidx = np.where(zi+zstep <= zCL)[0]
ABL_plumes = np.sort(np.array(RunList)[ABLidx])
penetrative_plumes = [i for i in RunList if i not in ABL_plumes]

#==============figure out distribution=======
prunCnt = len(PENidx)

zCLP = np.empty((prunCnt)) * np.nan
uBL = np.empty((prunCnt)) * np.nan
windRatio = np.empty((prunCnt)) * np.nan
zMax = np.empty((prunCnt)) * np.nan
zMaxGuess = np.empty((prunCnt)) * np.nan

profileModelled = np.empty((prunCnt,len(interpZ))) * np.nan
profileHalfModelled = np.empty((prunCnt,len(interpZ))) * np.nan

wF = ((g*Phi[PENidx]*(zCL[PENidx]-BLfrac*zi[PENidx]))/(thetaS[PENidx]*zi[PENidx]))**(1/3.)
wD = (g * zi[PENidx] * 0.13 / thetaS[PENidx])**(1/3.)

# exclude = ['W5F4R5TE','W5F4R6TE','W5F12R5TE','W5F13R5TE','W5F13R6TE']
exclude = []

for nCase,Case in enumerate(penetrative_plumes):
    if Case in exclude:
        continue
    else:
        zsidx = np.argmin(abs(interpZ - zS[PENidx][nCase]))
        ziidx = np.argmin(abs(interpZ - zi[PENidx][nCase]))
        zclidx = np.argmin(abs(interpZ - zCL[PENidx][nCase]))

        distribP = profile[PENidx][nCase,:]

        #create an interpolated profile of velocity
        if Case[-1:]=='T' or Case[-1:]=='E':
            levels = plume.lvltall
            if Case[-1:]=='E':
                pmlvl = np.arange(0,4001,40)        #extra tall domain for the concentrations only
            else:
                pmlvl=levels
        else:
            levels=plume.lvl

        U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

        interpU= interp1d(levels,U0,fill_value='extrapolate')
        U0interp = interpU(interpZ)
        uBL[nCase] = np.mean(U0interp[int(ziidx/2):ziidx])

        #find maximum
        cutoff = np.max(distribP) * 0.0015
        pmTopidx = np.argmin(abs(distribP[zclidx:] - cutoff)) + zclidx
        zMax[nCase] = interpZ[pmTopidx]
        sigmaTop = (zMax[nCase] - zCL[PENidx][nCase])/3.

        zMaxGuess[nCase] = 2*zCL[PENidx][nCase] - zS[PENidx][nCase]
        sigmaTopFair = (zMaxGuess[nCase] - zCL[PENidx][nCase])/3.

        #figure out the bottom half of the distribution
        windRatio[nCase] = uBL[nCase]/(wF[nCase] - wD[nCase])
        if wF[nCase]/wD[nCase] < 1.5:
            windFactor = uBL[nCase]/wF[nCase]
        else:
            windFactor = windRatio[nCase]

        if windFactor > 1:
            sigmaBottom = sigmaTop*windFactor
            sigmaBottomFair = sigmaTopFair * windFactor
        else:
            sigmaBottom = sigmaTop
            sigmaBottomFair = sigmaTopFair

        #make a fit with two gaussians
        y_pdf_Fair = norm.pdf(interpZ, zCL[PENidx][nCase], sigmaTopFair) # the normal pdf
        y_pdf_zMax = norm.pdf(interpZ, zCL[PENidx][nCase], sigmaTop) # the normal pdf
        y_pdf_Bottom_Fair = norm.pdf(interpZ, zCL[PENidx][nCase], sigmaBottomFair) # the normal pdf
        y_pdf_Bottom = norm.pdf(interpZ, zCL[PENidx][nCase], sigmaBottom) # the normal pdf

        profileModelled[nCase,zclidx:] = y_pdf_Fair[zclidx:]*np.max(distribP[zclidx:])/(np.max(y_pdf_Fair[zclidx:])*1000)
        profileHalfModelled[nCase,zclidx:] = y_pdf_zMax[zclidx:]*np.max(distribP[zclidx:])/(np.max(y_pdf_zMax[zclidx:])*1000)
        profileModelled[nCase,:zclidx+1] = y_pdf_Bottom_Fair[:zclidx+1]*np.max(distribP[:zclidx+1])/(np.max(y_pdf_Bottom_Fair[:zclidx+1])*1000)
        profileHalfModelled[nCase,:zclidx+1] = y_pdf_Bottom[:zclidx+1]*np.max(distribP[:zclidx+1])/(np.max(y_pdf_Bottom[:zclidx+1])*1000)

        plt.figure()
        plt.title('%s' %Case)
        plt.plot(distribP/1000,interpZ,label='smoke profile')
        ax = plt.gca()
        ax.set(xlabel='CWI concentration [ppm]',ylabel='height [m]',ylim=[0,3200])
        ax.fill_betweenx(interpZ, quartiles[PENidx][nCase,:,0]/1000,quartiles[PENidx][nCase,:,1]/1000, alpha=0.2,label='IQR')
        ax.axhline(y = interpZ[zclidx], ls='-.',c='C0',linewidth=1, label='z$_{CL}$')
        ax.axhline(y = zMax[nCase], ls='--', c='grey',label='$z_{top}$ LES' )
        ax.axhline(y = zMaxGuess[nCase],ls='--',c='C1',label='$z_{max}$ modelled' )
        plt.plot(profileHalfModelled[nCase,:],interpZ,c='grey',ls=':',label=r'based on LES z$_{top}$')
        plt.plot(profileModelled[nCase,:],interpZ,c='C1',label=r'based on modelled z$_{top}$')
        plt.legend()
        plt.savefig(plume.figdir + 'distribution/raw/pmProf%s.pdf' %Case)
        plt.close()

errorMax = zMax-zMaxGuess
topFit = linregress(zMaxGuess[np.isfinite(zMaxGuess)],zMax[np.isfinite(zMax)])
plt.figure()
plt.title('PREDICTING DISTRIBUTION TOP: R=%.2f' %topFit[2])
plt.scatter(zMaxGuess,zMax,c=zMax-zCL[PENidx])
plt.plot(np.arange(4500), np.arange(4500), c='grey',ls=':')
plt.gca().set(aspect='equal' , xlabel=r'$z_{CL} + z^\prime$ [m]', ylabel=r'$z_{top}$ from LES [m]',xlim=[700,4500],ylim=[700,4500])
plt.colorbar(label=r'$z^\prime$ [m]')
plt.tight_layout()
# plt.show()
plt.savefig(plume.figdir + 'distribution/TopPredictor.pdf' )
plt.close()

#error analysis
normprofile = np.empty((prunCnt,len(interpZ))) * np.nan
normprofileModelled = np.empty_like(normprofile) * np.nan
MAE = np.empty_like(zMax) * np.nan
MAE_BL = np.empty_like(zMax) * np.nan
MAE_FA = np.empty_like(zMax) * np.nan

normprofileHalfModelled = np.empty((prunCnt,len(interpZ))) * np.nan
MAE_half = np.empty_like(zMax) * np.nan
MAE_BL_half = np.empty_like(zMax) * np.nan
MAE_FA_half = np.empty_like(zMax) * np.nan

for nCase,Case in enumerate(penetrative_plumes):
    ziidx = np.argmin(abs(interpZ - zi[PENidx][nCase]))
    distribP = profile[PENidx][nCase,:]
    normtruth = distribP/np.nanmax(distribP)
    normmodel = profileModelled[nCase,:]/np.nanmax(profileModelled[nCase,:])
    normprofile[nCase,:] = normtruth
    normprofileModelled[nCase,:] = normmodel
    mae_subset = np.where((normtruth > 0.001) & (normmodel > 0.001))[0]
    MAE[nCase] = np.nanmean(abs(normtruth[mae_subset] - normmodel[mae_subset]))
    MAE_BL[nCase] = np.nanmean(abs(normtruth[mae_subset][mae_subset < ziidx] - normmodel[mae_subset][mae_subset < ziidx]))
    MAE_FA[nCase] = np.nanmean(abs(normtruth[mae_subset][mae_subset >= ziidx] - normmodel[mae_subset][mae_subset >= ziidx]))

    normhalfmodel = profileHalfModelled[nCase,:]/np.nanmax(profileHalfModelled[nCase,:])
    normprofileHalfModelled[nCase,:] = normhalfmodel
    mae_subset_half = np.where((normtruth > 0.001) & (normhalfmodel > 0.001))[0]
    MAE_half[nCase] = np.nanmean(abs(normtruth[mae_subset] - normmodel[mae_subset]))
    MAE_BL_half[nCase] = np.nanmean(abs(normtruth[mae_subset_half][mae_subset_half < ziidx] - normhalfmodel[mae_subset_half][mae_subset_half < ziidx]))
    MAE_FA_half[nCase] = np.nanmean(abs(normtruth[mae_subset_half][mae_subset_half >= ziidx] - normhalfmodel[mae_subset_half][mae_subset_half >= ziidx]))

plt.figure(figsize=(10,5))
plt.suptitle('NORMALIZED DISTRIBUTION MAE')
plt.subplot(121)
plt.title(r'LES-derived z$_{top}$')
plt.boxplot([MAE_BL_half[np.isfinite(MAE_BL_half)],MAE_FA_half[np.isfinite(MAE_FA_half)],MAE_half[np.isfinite(MAE_half)]], labels = ('ABL','FREE ATM','TOTAL'))
plt.gca().set(ylabel='MAE (normalized concentration)',ylim=[0,0.4])
plt.subplot(122)
plt.title(r'Modelled z$_{top}$')
plt.boxplot([MAE_BL[np.isfinite(MAE_BL)],MAE_FA[np.isfinite(MAE_FA)],MAE[np.isfinite(MAE)]], labels = ('ABL','FREE ATM','TOTAL'))
plt.gca().set(ylabel='MAE (normalized concentration)',ylim=[0,0.4])
plt.tight_layout()
plt.savefig(plume.figdir + 'distribution/MAEdistribution.pdf')
plt.close()

plt.figure()
plt.suptitle('ABL MAE')
plt.scatter(plume.read_tag('W', RunList)[PENidx],MAE_BL,c=uBL/(wF-wD))
plt.colorbar(label=r'$z_i$ [m]')
plt.gca().set(xlabel='wind [m/s]',ylabel='ABL MAE (normalized concentration)')
plt.savefig(plume.figdir + 'distribution/ABL_MAEvsWind.pdf')
plt.close()
