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

#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

testPortion = 0.2       #portion of data to reserve for testing the model
trials = 10             #number of times to rerun the model
#=================end of input===============

# RunList = [i for i in plume.tag if i not in plume.exclude_runs]
RunList = ['W5F1R0']

runCnt = len(RunList)
g = 9.81
si = 3

#saved variables
r = np.empty((runCnt)) * np.nan                 #fireline depth
Ua = np.empty((runCnt)) * np.nan                #ambient wind
zi = np.empty((runCnt)) * np.nan                #BL height
zCL = np.empty((runCnt)) * np.nan               #centerline height
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature
Phi = np.empty((runCnt)) * np.nan               #cumulative fire heat
FI = np.empty((runCnt)) * np.nan                #total 2D fire heating
BL = np.empty((runCnt,len(plume.lvl)-1))        #storage for BL
FlaggedCases = []                               #for storage of indecies of anomalous runs


#======================repeat main analysis for all runs first===================

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
    smoothCenterline = savgol_filter(centerline, 31, 3) # window size 31, polynomial order 3

    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 101, 3) # window size 101, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.05 and nX > np.nanargmax(smoothPM) else False for nX in range(dimX-1) ]
    stablePM = pm[:,1:][:,stablePMmask]
    stableProfile = np.mean(stablePM,1)
    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)

    #define source (r and H)------------------------
    #raduis using full 2D average -
    masked_flux = ma.masked_less_equal(csdict['ghfx2D'], 0)
    cs_flux = np.nanmean(masked_flux,1)                         #get cross section for each timestep
    fire = []
    xmax = np.argmax(cs_flux,axis=1)                            #get max for each timestep
    for nP, pt in enumerate(xmax[plume.ign_over:]):             #excludes steps containing ignition
        subset = cs_flux[plume.ign_over+nP,pt-plume.wi:pt+plume.wf]     #set averaging window around a maximum
        fire.append(subset)
    burning = np.sum(np.sum(csdict['ghfx2D'][plume.ign_over:,:,:],1),1) #get total heat flux from entire fire area

    meanFire = np.nanmean(fire,0)
    ignited = np.array([i for i in meanFire if i > 0.5])
    r[nCase] = len(ignited) * plume.dx
    Phi[nCase] = np.trapz(ignited, dx = plume.dx) * 1000 / ( 1.2 * 1005)
    FI[nCase] = np.mean(burning)*  1000 / ( 1.2 * 1005)

    #injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])
    zCLidx = int(np.mean(ctrZidx[1:][stablePMmask]))
    dT = T0[1:]-T0[0:-1]
    Omega[nCase] = np.trapz(dT[si+1:zCLidx], dx = plume.dz)
    BL[nCase,:] = dT

    #highlight weird plumes that don't reach top of boundary layer
    if Omega[nCase] < 0 :
        print('\033[93m' + 'Omega: %0.2f ' %Omega[nCase] + '\033[0m')
        print('\033[93m' + 'Hard overwrite (see VelScale_InjHeight.py): Omega = Omega[zi]' + '\033[0m' )
        ziIdx = np.where(plume.lvl==zi[nCase])[0][0]
        Omega[nCase] = np.trapz(dT[si+1:ziIdx], dx = plume.dz)
        FlaggedCases.append(nCase)
    Ua[nCase] = np.mean(U0[si:zCLidx])


#======================train and test regression model===================

#linear regression using all data
wStar = (g*Phi*zi/(Omega))**(1/3.)

slopeALL, interceptALL, r_valueALL, p_valueALL, std_errALL = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])

Rstore = np.empty((trials)) * np.nan
ModelError = []
TrueTrialZcl = []

for nTrial in range(trials):
    #split runs into train and test datasets
    TestFlag = np.random.binomial(1,0.2,runCnt)
    testCnt = sum(TestFlag)

    #linear regression using training data
    slope, intercept, r_value, p_value, std_err = linregress(wStar[TestFlag==0][np.isfinite(wStar[TestFlag==0])],zCL[TestFlag==0][np.isfinite(wStar[TestFlag==0])])
    print('Sum of residuals using TRAINING data: %0.2f' %r_value)
    Rstore[nTrial] = r_value

    fig = plt.figure()
    plt.suptitle('REGRESSION MODEL: ALL [R=%0.2f] vs TRAIN DATA [R=%0.2f]' %(r_valueALL, r_value))
    ax=plt.gca()
    plt.scatter(wStar[TestFlag==0], zCL[TestFlag==0], c='C2', label='training data')
    plt.scatter(wStar[TestFlag==1], zCL[TestFlag==1], c='C1', label='test data')
    plt.plot(wStar, interceptALL + slopeALL*wStar, c='grey', label='all data')
    plt.plot(wStar, intercept + slope*wStar, c='C2', label='training data regression fit')
    ax.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
    plt.legend()
    plt.savefig(plume.figdir + 'injectionModel/trials/ALLvsTRAINdataTrial%s.pdf' %nTrial )
    plt.show()
    plt.close()

    '''
    Solve a system of equations:
    (1) zCL = m*wStar + b
    (2) wStar = [g * Phi * zi / Omega]**1/3
    '''

    #solve numerically
    zCLmodel = np.empty((testCnt)) * np.nan

    for nTest in range(testCnt):
        toSolve = lambda z : z - intercept - slope * (g*Phi[TestFlag==1][nTest]*zi[TestFlag==1][nTest]/(np.trapz(BL[TestFlag==1][nTest][si+1:int(z/plume.dz)], dx = plume.dz)))**(1/3.)
        z_initial_guess = zi[TestFlag==1][nTest]                    #make initial guess BL height
        z_solution = fsolve(toSolve, z_initial_guess)

        zCLmodel[nTest] = z_solution
        print('%s solution is zCL = %0.2f' % (np.array(RunList)[TestFlag==1][nTest],z_solution))
        print('...True value: %0.2f ' %zCL[TestFlag==1][nTest])
    error = zCLmodel -  zCL[TestFlag==1]
    ModelError.append(error)
    TrueTrialZcl.append(zCL[TestFlag==1])

print('Sum of residuals using ALL data: %0.2f' %r_valueALL)
print('\033[93m' + 'Linear model equation using ALL data: zCL = %.3f Wf* + %.3f '  %(slopeALL,interceptALL)+ '\033[0m' )

#======================plot model stability===================

flatTrueTrialZcl  = np.concatenate(TrueTrialZcl)                #flatten test array
flatModelError = np.concatenate(ModelError)

plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
ax0 = plt.subplot(gs[0])
plt.title('TRIAL ERROR DISTRIBUTIONS')
plt.boxplot(ModelError)
plt.hlines(0,0,10,colors='grey',linestyles='dashed')
ax0.set(xlabel='trial no.', ylabel='error in zCL [m]')
ax1 = plt.subplot(gs[1])
plt.title('TOTAL ERROR DISTRIBUTION')
plt.boxplot(flatModelError)
plt.hlines(0,0,10,colors='grey',linestyles='dashed')
ax1.set(xlabel='all runs', ylabel='error in zCL [m]')
ax2 = plt.subplot(gs[2])
plt.title('PREDICTION INTERVAL')
m, b, r, p, std = linregress(flatTrueTrialZcl,flatModelError)
plt.scatter(flatTrueTrialZcl,flatModelError)
plt.plot(flatTrueTrialZcl, m*flatTrueTrialZcl + b, color='grey')
ax2.set(xlabel='plume height [m]', ylabel='model error [m]')
ax3 = plt.subplot(gs[3])
plt.title('R-VALUE SENSITIVITY')
plt.hist(Rstore, bins=5)
ax3.set(xlabel='R-value',ylabel='count' )
plt.tight_layout()
plt.savefig(plume.figdir + 'injectionModel/ModelSensitivity.pdf')
plt.show()
plt.close()
