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


#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time

testPortion = 0.2       #portion of data to reserve for testing the model
trials = 10             #number of times to rerun the model

g = 9.81                #gravity constant
zs= 500                 #surface layer height in m
zstep = 20              #height interpolation step

#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]         #load a list of cases
runCnt = len(RunList)                           #count number of cases

#set up interpolated vertical profile with 5m vertical step
interpZ = np.arange(0, plume.lvltall[-1], zstep)
si = int(zs/zstep)


#storage for variables
zi = np.empty((runCnt)) * np.nan                #BL height (m)
zCL = np.empty((runCnt)) * np.nan               #smoke injection height (m)
Phi = np.empty((runCnt)) * np.nan               #cumulative fire heat  (K m^2/s)
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature (Km) - the questionable denominator term
Omega2 = np.empty((runCnt)) * np.nan            #ALTERNATIVE TERM: simply potential temperature change between surface layer and injection height (K)
varTest = np.empty((runCnt)) * np.nan
sounding = np.empty((runCnt,len(interpZ))) * np.nan         #storage for interpolated soundings
gradT0interp = np.empty((runCnt,len(interpZ)-1)) * np.nan   #storage for temperature gradient
zCLerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
fireEnd = np.empty((runCnt)) * np.nan
zSTD = np.empty((runCnt)) * np.nan
Wcum = np.empty((runCnt)) * np.nan

FlaggedCases = []                               #for storage of indecies of anomalous runs

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
    if Case[-2:]=='6T':
        interpT = interp1d(plume.lvltall, T0,fill_value='extrapolate')
    else:
        interpT= interp1d(plume.lvl, T0,fill_value='extrapolate')
    T0interp = interpT(interpZ)


    #mask plume with cutoff value---------------------------------
    dimT, dimZ, dimX = np.shape(csdict['temp'])     #get shape of data
    zi[nCase] = plume.get_zi(T0)                    #calculate BL height
    if Case[-2:]=='6T':
        zi[nCase] = 2211
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= plume.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells


    #locate centerline
    ctrZidx = pm.argmax(0)                          #locate maxima along height
    ctrXidx = pm.argmax(1)                          #locate maxima downwind
    pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])    #get concentration along the centerline

    xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)           #get location of maximum centerline height
    centerline = ma.masked_where(plume.lvltall[ctrZidx] == 0, plume.lvltall[ctrZidx])               #make sure centerline is only calculated inside the plume
    smoothCenterline = savgol_filter(centerline, 71, 3)             # smooth centerline height (window size 31, polynomial order 3)


    #calculate concentration changes along the centerline
    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, 101, 3) # window size 101, polynomial order 3
    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                            abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                            nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                            nX > np.nanargmax(smoothPM) and\
                            nX > np.nanargmax(centerline) +10 and\
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
    stableProfile = np.mean(stablePM,1)

    pmQ1 = np.percentile(stablePM,25,axis = 1)
    pmQ3 = np.percentile(stablePM,75,axis = 1)

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

    # #comparison with last frame (to make sure averaging is a good idea)
    # plt.figure()
    # plt.title(RunList[nCase])
    # plt.plot(meanFire, label='mean')
    # plt.plot(subset,label='last')
    # plt.legend()
    # plt.savefig(plume.figdir + 'temp/fireline%s.pdf' %Case)
    # plt.close()

    Phi[nCase] = np.trapz(ignited, dx = plume.dx) * 1000 / ( 1.2 * 1005)    #calculate Phi by integrating kinematic heat flux along x (Km2/s)

    #calculate injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
    dT = (T0[1:]-T0[0:-1])/plume.dz                                        #calculate potential temperature change (K)
    zSTD[nCase] = np.std(smoothCenterline[1:][stablePMmask])

    siBL = np.argmin(abs(interpZ - zi[nCase]*2/3))
    sounding[nCase,:] = T0interp
    dTinterp = (T0interp[1:] - T0interp[0:-1])/zstep
    gradT0interp[nCase,:] = dTinterp
    # Omega[nCase] = np.trapz(dTinterp[si:zCLidx], dx = zstep)     #calculate the denominator term by integrating temperature change with height, excluding surface layer (Km)
    Omega[nCase] = np.trapz(dTinterp[siBL:zCLidx], dx = zstep)     #calculate the denominator term by integrating temperature change with height, excluding surface layer (Km)

    w = csdict['w'][-1,:,:]
    Wcum[nCase] = np.sum(np.max(w,1)[:ymax+5])
    # if zCL[nCase]<zi[nCase]:
    #     print('INJECTION IS BELOW THE BL HEIGHT!!!!!!!!!!!!!!!')
    #     Omega[nCase] = np.trapz(dTinterp[si:np.argmin(abs(interpZ - zi[nCase]))], dx = plume.zstep)


    #alternative Omega calculation (simply temperature change)
    # Omega2[nCase] = (T0interp[zCLidx] - T0interp[si])                      # change in potential temperature between surface layer and injection layer (K)
    Omega2[nCase] = (T0interp[zCLidx] - np.mean(T0interp[siBL-2:siBL+2]))                      # change in potential temperature between surface layer and injection layer (K)

    zCLidxNONINTERP = np.argmin(abs(plume.lvltall - zCL[nCase]))

    varTest[nCase] = (fxmax[-1]+1)*plume.dx - 1000
    # #highlight weird plumes that don't reach top of boundary layer
    # if Omega[nCase] < 0 :
    #     print('\033[93m' + 'Omega: %0.2f ' %Omega[nCase] + '\033[0m')
    #     print('\033[93m' + 'Hard overwrite (see VelScale_InjHeight.py): Omega = Omega[zi]' + '\033[0m' )
    #     ziIdx = np.where(plume.lvl==zi[nCase])[0][0]
    #     Omega[nCase] = np.trapz(dT[si+1:ziIdx], dx = plume.dz)
    #     FlaggedCases.append(nCase)
#======================compare model formulations========================
#define wf* (as per original 'wrong' formulation)
# wStar = (g*Phi*(zi-zs)/(Omega))**(1/3.)

wStar = (g*Phi*(zi-zs)/(Omega))**(1/3.)

ddorf = (g*0.2*zi/sounding[:,25])**(1/3.)

plt.scatter(wStar,Wcum)
plt.plot(np.arange(1000),np.arange(1000))
plt.show()

#do linear regression using all data
slopeALL, interceptALL, r_valueALL, p_valueALL, std_errALL = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])

#using formulation suggested by roland
wStar2 = (g * Phi * (zi-zs)/ (Omega2))**(1/3.)

slope2, intercept2, r_value2, p_value2, std_err2 = linregress(wStar2[np.isfinite(wStar2)],zCL[np.isfinite(wStar2)])


#make scatterplot comparisons
plt.figure(figsize=(12,5))
ax1 = plt.subplot(121)
plt.title('MATHEMATICAL FORMULATION: [R=%0.2f]' %r_value2)
plt.scatter(wStar2, zCL, c=plume.read_tag('W',RunList), label=r'$w_{f*} = \frac{g \cdot (z_i-z_s) \cdot \Phi \cdot \epsilon}{(\theta_{z_{CL}} - \theta_{z_s})}$')
ax1.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
plt.colorbar().set_label('ambient wind [m/s]')
plt.legend(fontsize=12)
ax1 = plt.subplot(122)
plt.title('NUMERICAL FORMULATION: [R=%0.2f]' %r_valueALL)
plt.scatter(wStar, zCL, c=Phi, cmap = plt.cm.Spectral_r, label=r'$w_{f*} = \frac{g \cdot (z_i - z_s) \cdot \Phi \cdot \epsilon }{\int_{z_s}^{z_{cl}} d\theta  dz}}$')
plt.plot(wStar, interceptALL + slopeALL*wStar,c='grey')
plt.colorbar().set_label('$\Phi$ [Km$^2$/s]')
plt.legend(fontsize=12)
ax1.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
plt.tight_layout()
plt.savefig(plume.figdir + 'injectionModel/CompareFormulations.pdf')
plt.show()
plt.close()

#======================train and test regression model===================


#create storage arrays for R values, modelled zCL, model error and trail subsets of true zCL derived from data
Rstore = np.empty((trials)) * np.nan
ModelError = []
TrueTrialZcl = []
TrialFuel = []
TrialVar = []
TrialName = []

for nTrial in range(trials):
    #split runs into train and test datasets
    TestFlag = np.random.binomial(1,testPortion,runCnt)                 #pick a random set of approximately 20% of data
    testCnt = sum(TestFlag)         #count how many runs ended up as test dataset

    #linear regression using training data subset only
    slope, intercept, r_value, p_value, std_err = linregress(wStar[TestFlag==0][np.isfinite(wStar[TestFlag==0])],zCL[TestFlag==0][np.isfinite(wStar[TestFlag==0])])
    print('Sum of residuals using TRAINING data: %0.2f' %r_value)
    Rstore[nTrial] = r_value        #store trial value

    #plot individual trial results
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
    NOW APPLY THE MODEL
    Solve a system of equations:
    (1) zCL = m*wStar + b
    (2) wStar = [g * Phi * zi / Omega]**1/3
    Use values m and b calcuated for the training dataset only and apply them to the test dataset
    '''

    #solve numerically for each run belonging to the Test subgroup
    zCLmodel = np.empty((testCnt)) * np.nan

    for nTest in range(testCnt):
        testIdx = np.where(TestFlag==1)[0][nTest]                   #get index of test run in the LES subset
        #substituting equation (2) into (1) above
        toSolve = lambda z : z - intercept - slope * \
                            ( (g*Phi[testIdx]*(zi[testIdx] - zs))/ \
                            (np.trapz(gradT0interp[testIdx][si:int(z/zstep)], dx=zstep) ))**(1/3.)

        z_initial_guess = zi[testIdx]                    #make initial guess BL height
        z_solution = fsolve(toSolve, z_initial_guess)               #solve

        zCLmodel[nTest] = z_solution                                #store the solution
        print('%s solution is zCL = %0.2f' % (np.array(RunList)[testIdx],z_solution))
        print('...True value: %0.2f ' %zCL[testIdx])
    error = zCLmodel -  zCL[TestFlag==1]                            #calculate error between model and 'truth'
    ModelError.append(error)                                        #store model error
    TrueTrialZcl.append(zCL[TestFlag==1])                           #store true subset
    category = plume.read_tag('F',np.array(RunList)[TestFlag==1])
    TrialFuel.append(category)
    TrialVar.append(zCL[TestFlag==1]-zs)
    TrialName.append(np.array(RunList)[TestFlag==1])

print('Sum of residuals using ALL data: %0.2f' %r_valueALL)
print('\033[93m' + 'Linear model equation using ALL data: zCL = %.3f Wf* + %.3f '  %(slopeALL,interceptALL)+ '\033[0m' )

#======================plot model stability===================

flatTrueTrialZcl  = np.concatenate(TrueTrialZcl)                #flatten test array of injection heights
flatModelError = np.concatenate(ModelError)                     #flatten model error
flatTrialFuel = np.concatenate(TrialFuel)
flatTrialVar = np.concatenate(TrialVar)
flatTrialName = np.concatenate(TrialName)


#plot model sensitivity
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
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel='all runs', ylabel='error in zCL [m]')
ax2 = plt.subplot(gs[2])
plt.title('PREDICTION INTERVAL')
m, b, r, p, std = linregress(flatTrueTrialZcl,flatModelError)
plt.scatter(flatTrueTrialZcl,flatModelError)
plt.plot(flatTrueTrialZcl, m*flatTrueTrialZcl + b, color='grey')
plt.hlines(0,500,2000,colors='grey',linestyles='dashed')
ax2.set(xlabel='plume height [m]', ylabel='model error [m]',xlim=[500,2000])
ax3 = plt.subplot(gs[3])
plt.title('R-VALUE SENSITIVITY')
plt.hist(Rstore, bins=5)
ax3.set(xlabel='R-value',ylabel='count' )
plt.tight_layout()
plt.savefig(plume.figdir + 'injectionModel/ModelSensitivity.pdf')
plt.show()
plt.close()


plt.figure()
plt.title('ERROR AS A FUNCTION OF FUEL (TRIALS)')
plt.scatter(flatTrialFuel,flatModelError,c=flatTrialVar)
plt.hlines(0,0,14,colors='grey',linestyles='dashed')
ax = plt.gca()
# for i, txt in enumerate(flatTrialName):
#     ax.annotate(txt, (flatTrialFuel[i], flatModelError[i]),fontsize=6)
ax.set(xlabel='fuel category', ylabel='error [m]',ylim=[-100,150])
plt.colorbar().set_label('$z_{CL} - z_s$ [m]')
plt.savefig(plume.figdir + 'injectionModel/FuelvsErrorHeight_TRIALS.pdf')
plt.show()
plt.close()


#============================model sensitivity and entrainment============

#define wf* (as per original 'wrong' formulation)
wStar = (g*Phi* (zi-zs)/(Omega))**(1/3.)
# wStar = (g*Phi* (zi-zs)/(Omega ))**(1/3.)

#do linear regression using all data
slopeALL, interceptALL, r_valueALL, p_valueALL, std_errALL = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])
print(r_valueALL)
plt.scatter(wStar,zCL,c=zSTD, cmap=plt.cm.rainbow )
plt.colorbar()
ax = plt.gca()
for i, txt in enumerate(RunList):
    ax.annotate(txt, (wStar[i], zCL[i]),fontsize=6)
plt.savefig(plume.figdir + 'injectionModel/InjectionModelSTD.pdf')

plt.show()



for nCase,Case in enumerate(RunList):
    toSolveCase = lambda z : z - interceptALL - slopeALL * \
                            ( (g*Phi[nCase]*(zi[nCase] - zs))/ \
                            (np.trapz(gradT0interp[nCase][si:int(z/zstep)], dx=zstep) ))**(1/3.)

    z_initial_guess = zi[nCase]                    #make initial guess BL height
    z_solution = fsolve(toSolveCase, z_initial_guess)               #solve
    zCLerror[nCase] = z_solution - zCL[nCase]                                #store the solution

plt.figure()
plt.title('ERROR AS A FUNCTION OF zCL (ALL)')
plt.scatter(zCL,zCLerror,c=varTest)
# mBIAS, bBIAS, r, p, std = linregress(zCL,zCLerror)
ax = plt.gca()
for i, txt in enumerate(RunList):
    ax.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
plt.plot(zCL, m*zCL + b, color='grey')
plt.show()
# plt.close()



plt.figure()
plt.title('ERROR AS A FUNCTION OF FUEL (ALL)')
plt.scatter(plume.read_tag('F',RunList),zCLerror,c=zCL)
plt.hlines(0,0,14,colors='grey',linestyles='dashed')
ax.set(xlabel='fuel category', ylabel='error [m]')
plt.colorbar().set_label('$z_{CL}$ [m]')
plt.savefig(plume.figdir + 'injectionModel/FuelvsErrorHeight_ALL.pdf')
plt.show()
# plt.close()

plt.figure()
plt.title('PRE-IGNITION ATMOSPHERIC PROFILES')
leg_handles = []
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)
for R in set(Rtag):
    for Case in sounding[Rtag==R]:
        lR = plt.plot(Case, interpZ, color='C%s' %R, linewidth=1, label='R%s' %R)
    leg_handles.extend(lR)
plt.gca().set(xlabel='potential temperature [K]',ylabel='height [m]',xlim=[280,330],ylim=[0,2800])
plt.legend(handles=leg_handles)
plt.savefig(plume.figdir + 'T0profiles.pdf')
plt.show()
