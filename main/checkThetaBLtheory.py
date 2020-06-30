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
Phi2 = np.empty((runCnt)) * np.nan
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature (Km) - the questionable denominator term
Omega2 = np.empty((runCnt)) * np.nan


sounding = np.empty((runCnt,len(interpZ))) * np.nan         #storage for interpolated soundings
gradT0interp = np.empty((runCnt,len(interpZ)-1)) * np.nan   #storage for temperature gradient
wM = np.empty((runCnt)) * np.nan
wCum = np.empty((runCnt)) * np.nan
wCum2 = np.empty((runCnt)) * np.nan
wCum3 = np.empty((runCnt)) * np.nan
wCum4 = np.empty((runCnt)) * np.nan
wCum5 = np.empty((runCnt)) * np.nan

depth = np.empty((runCnt)) * np.nan

thetaCL = np.empty((runCnt)) * np.nan
thetaS = np.empty((runCnt)) * np.nan
thetaZI = np.empty((runCnt)) * np.nan

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
        interpT = interp1d(plume.lvltall, T0,fill_value='extrapolate')
    else:
        interpT= interp1d(plume.lvl, T0,fill_value='extrapolate')
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
    filter_window = max(int(plume.read_tag('W',[Case])*10+1),51)
    smoothCenterline = savgol_filter(centerline, filter_window, 3)             # smooth centerline height (window size 31, polynomial order 3)


    #calculate concentration changes along the centerline
    dPMdX = pmCtr[1:]-pmCtr[0:-1]
    smoothPM = savgol_filter(dPMdX, filter_window, 3) # window size 101, polynomial order 3

    stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                            abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                            nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                            nX > np.nanargmax(smoothPM) and\
                            nX > np.nanargmax(centerline) +10 and
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
    Phi[nCase] = np.trapz(ignited, dx = plume.dx) * 1000 / ( 1.2 * 1005)    #calculate Phi by integrating kinematic heat flux along x (Km2/s)
    Phi2[nCase] = np.sum(ignited) * 1000 / ( 1.2 * 1005)
    depth[nCase] = len(ignited) * plume.dx
    #calculate injection height variables ---------------------------
    zCL[nCase] = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
    zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
    ziidx = np.argmin(abs(interpZ - zi[nCase]))
    dT = (T0[1:]-T0[0:-1])/plume.dz                                        #calculate potential temperature change (K)

    siBL = np.argmin(abs(interpZ - zi[nCase]*2/3))
    meanBLtemp,siBLtemp = np.mean(T0interp[10:ziidx]), T0interp[siBL]
    print(meanBLtemp - siBLtemp)
    sounding[nCase,:] = T0interp
    dTinterp = (T0interp[1:] - T0interp[0:-1])/zstep
    gradT0interp[nCase,:] = dTinterp
    Omega[nCase] = np.trapz(dTinterp[si:zCLidx], dx = zstep)     #calculate the denominator term by integrating temperature change with height, excluding surface layer (Km)
    Omega2[nCase] = np.trapz(dTinterp[siBL:zCLidx],dx=zstep)


    w = csdict['w'][-1,:,:]
    Wctr = np.array([w[nZ, ctrXidx[nZ]] for nZ in range(dimZ)])    #get concentration along the centerline
    Wmax = np.max(w,1)
    # Wcum[nCase] = np.sum(np.max(w,1)[:ymax+5])
    if Case[-1:]=='T' or Case[-1:]=='E':
        interpWctr = interp1d(plume.lvltall, Wctr,fill_value='extrapolate')
        interpWmax = interp1d(plume.lvltall, Wmax,fill_value='extrapolate')
    else:
        interpWctr= interp1d(plume.lvl, Wctr,fill_value='extrapolate')
        interpWmax= interp1d(plume.lvl, Wmax,fill_value='extrapolate')
    Wctrinterp = interpWctr(interpZ)
    Wmaxinterp = interpWmax(interpZ)
    wM[nCase] = max(Wmax)
    wCum[nCase] = np.sum(Wmaxinterp[25:zCLidx])
    wCum2[nCase] = np.sum(Wctrinterp[25:zCLidx])
    wCum3[nCase] = np.trapz(Wmaxinterp[25:zCLidx], dx=40.)
    wCum4[nCase] = np.trapz(Wctrinterp[25:zCLidx], dx=40.)

    #do theta testing
    temperature = csdict['temp'][-1,:,:]
    Tctr = np.array([temperature[nZ, ctrXidx[nZ]] for nZ in range(dimZ)])    #get concentration along the centerline
    Tmax = np.max(temperature,1)
    if Case[-1:]=='T' or Case[-1:]=='E':
        interpTctr = interp1d(plume.lvltall, Tctr,fill_value='extrapolate')
        interpTmax = interp1d(plume.lvltall, Tmax,fill_value='extrapolate')
    else:
        interpTctr= interp1d(plume.lvl, Tctr,fill_value='extrapolate')
        interpTmax= interp1d(plume.lvl, Tmax,fill_value='extrapolate')
    Tctrinterp = interpTctr(interpZ)
    Tmaxinterp = interpTctr(interpZ)
    thetaCL[nCase] = sounding[nCase,zCLidx]
    thetaS[nCase] = sounding[nCase,25]
    thetaZI[nCase] = sounding[nCase,siBL]

    #vertical concentration slice at donwind locations of wmax and qmax
    plt.figure(figsize=(10,4))
    plt.suptitle('%s' %Case)
    plt.subplot(121)
    ax1 = plt.gca()
    plt.title('PROFILES OF VERTICAL VELOCITY')
    plt.plot(Wctrinterp,interpZ,'.-',label='$w_{PMmax}$')
    plt.plot(Wmaxinterp,interpZ,'k.-',label='$w_{max}$')
    plt.axhline(y = zi[nCase], ls=':', c='darkgrey', label='zi')
    plt.axhline(y = zCL[nCase],ls='--', c='red',label='z$_{CL}$')
    ax1.set(xlabel = 'velocity [m/s]', ylabel='height [m]',ylim = [0,3200] )
    plt.legend()

    plt.subplot(122)
    plt.title('PLUME vs AMBIENT TEMPERATURE')
    ax2 = plt.gca()
    plt.plot(sounding[nCase,:], interpZ, label='pre-ignition profile',c='lightblue')
    plt.plot(Tctrinterp,interpZ,c = 'orange',label='in-plume T$_{PMmax}$',alpha = 0.5)
    plt.plot(Tmaxinterp,interpZ,c = 'maroon',label='in-plume T$_{max}$')
    plt.axhline(y = zi[nCase], ls=':', c='darkgrey', label='zi')
    plt.axhline(y = zCL[nCase],ls='--', c='red',label='z$_{CL}$')
    ax2.set(xlabel = 'temperature [K]', ylabel='height [m]' ,xlim = [285,330],ylim = [0,3200])
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig(plume.figdir + 'profiles/profiles_%s.pdf' %Case)
    plt.close()


#======================compare model formulations========================
tau500 = 1/ np.sqrt(g*Omega/(thetaS * (zCL-500)))
A500 = np.log(thetaCL/thetaS)**(1/3.)

tauZI = 1/ np.sqrt(g*Omega2/(thetaZI * (zCL-zi*(2/3.))))
AZI = np.log(thetaCL/thetaZI)**(1/3.)
#define wf* (as per original 'wrong' formulation)
wStar = (1/40.)**(1/3.) *(g*Phi2*(zCL-500)*(3/2.)/(Omega))**(1/3.)              #our standard definition as sum(H)
wStar2 =  A500* (1/40.)**(1/3.) * (g*Phi2*(zCL-500)*(3/2.)/(Omega))**(1/3.)     #proper form (eqn 20)
wStar3 = (1/40.)**(1/3.) *(g*Phi2*(zCL-500)*(3/2.)/(thetaS))**(1/3.)            #proper form (eqn 18)

wStar4 = tauZI*(g*Phi*(zCL-zi*(2/3.))*(3/2.)/(thetaZI * zi))**(1/3.)       #attempt at zi_dependent threshold

wStar5 = A500*tau500*(1/40.)*((g*Phi2*(zCL-500))*(3/2.)/(Omega))**(1/3.)        #eqn 18 * 1/N
wStar6 = tau500*((1/40.)**(1/3.))*(g*Phi2*(zCL-500)*(3/2.)/(thetaS))**(1/3.)    #eqn 20 * 1/N
wStar7 = tau500*(g*Phi*(zCL-500)*(3/2.)/(thetaS * zi))**(1/3.)    #eqn 20 * 1/N


#do linear regression using all data
wStarFit = linregress(wStar,zCL)
wStarFit2 = linregress(wStar2,zCL)
wStarFit3 = linregress(wStar3,zCL)
wStarFit4 = linregress(wStar4+zi*(2/3),zCL)
wStarFit5 = linregress(wStar5,zCL)
wStarFit6 = linregress(wStar6,zCL)
wStarFit7 = linregress(wStar7,zCL)


wStar8 =  (3/4) *tauZI*(g*Phi*(zCL-zi*(2/3.))*(3/2.)/(thetaZI*zi))**(1/3.)       #attempt at zi_dependent threshold
wStarFit8 = linregress(wStar8+zi*(2/3),zCL)
print(wStarFit8)




print(wStarFit)
print(wStarFit2)
print(wStarFit3)
print(wStarFit4)
print(wStarFit5)
print(wStarFit6)
print(wStarFit7)
print(wStarFit8)

# plt.figure()
# plt.scatter(wStar6, zCL,c=plume.read_tag('W',RunList))
# plt.plot(wStar6,wStar6+500)
# plt.show()
#
# plt.figure()
# plt.scatter(wStar7, zCL,c=plume.read_tag('W',RunList))
# plt.plot(wStar7,wStar7+500)
# plt.show()

plt.figure()
plt.title('MODELLED SMOKE INJECTION HEIGHTS')
ax = plt.gca()
plt.scatter(wStar8+zi*(2/3),zCL,c=Phi,cmap =plt.cm.plasma)
ax.set(ylabel = r'$z_{CL}$ [m]', xlabel = r'$\frac{3}{4}\tau_* w_{f*} + \frac{2}{3}z_i$ [m]',xlim = [400,3200], ylim = [400,3200])
# for i, txt in enumerate(RunList):
#     ax.annotate(txt, (wStar4[i]+zi[i]*(2/3), zCL[i]),fontsize=6)
plt.colorbar(label=r'fireline intensity [K m$^2$/s]')
plt.plot(np.sort(wStar8+zi*(2/3)),wStarFit8[0]* np.sort(wStar8+zi*(2/3)) + wStarFit8[1], color='black', label='linear regression fit')
plt.plot(np.sort(wStar8+zi*(2/3)),np.sort(wStar8+zi*(2/3)), linestyle = 'dashed', color='grey', label='unity line')
plt.legend()
plt.savefig(plume.figdir + 'injectionModel/NewInjectionTheory.pdf')

plt.show()

# plt.scatter(wStar, zCL)
# plt.show()
#
#
# plt.figure()
# # plt.scatter(wStar,zCL)
# plt.scatter(wStar2,zCL,c = plume.read_tag('W',RunList)/ depth)
# plt.plot(wStar2,wStarFit2[0]*wStar2 + wStarFit2[1])
# # plt.gca().set(aspect='equal')
# plt.show()
#
# plt.scatter(wStar,(wCum2*depth)/(zCL-500))
# plt.gca().set(aspect='equal')
# plt.show()



# Omega500 = np.copy(Omega)
# plt.figure()
# ax1 = plt.gca()
# ax2 = plt.twinx()
# for Zr in np.arange(100,501,100):
#     wStartest = (((3/2.)*g*Phi* (zCL-Zr))/(vars()['Omega%s' %Zr]))**(1/3.)
#     c1, Zr_out, r,p,std = linregress(wStartest,zCL)
#     print(r)
#     ax1.scatter(Zr,Zr_out,c='C1')
#     ax2.scatter(Zr,r,c='C2')
# l1 = ax1.scatter(Zr,Zr_out,c='C1', label='Zr from fit')
# l2 = ax2.scatter(Zr,r,c='C2', label='R of the fit')
# ax1.set(xlabel='Zr_in [m]',aspect='equal')
# ax1.set_ylabel('Zr_out [m]', color='C1')
# ax2.set_ylabel('R value',color='C2')
# plt.title('Zi CUTOFF: %s' %500)
# plt.legend(handles=[l1,l2])
# plt.show()


# plt.figure()
# plt.scatter(wStar,wCum)
# plt.gca().set(aspect='equal')
# plt.show()

# plt.scatter(wStar,Wcum)
# plt.plot(np.arange(1000),np.arange(1000))
# plt.show()
# plt.close()



# #make scatterplot comparisons
# plt.figure(figsize=(12,5))
# plt.subplot(121)
# ax1 = plt.gca()
# plt.title('Wf* using zi: [R=%0.2f]' %fitZi[2])
# plt.scatter(wStar_zi, zCL, c=Phi, cmap = plt.cm.Spectral_r, label=r'$w_{f*} = \frac{g \cdot (z_i - z_s) \cdot \Phi \cdot \epsilon }{(\theta_{CL} - \theta_{BL})}}}$')
# plt.plot(wStar_zi, fitZi[1] + fitZi[0]*wStar_zi,c='grey')
# plt.colorbar().set_label('$\Phi$ [Km$^2$/s]')
# plt.legend(fontsize=12)
# ax1.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
# plt.subplot(122)
# ax2 = plt.gca()
# plt.title('Wf* using zCL: [R=%0.2f]' %fitZi[2])
# plt.scatter(wStar_zCL, zCL, c=Phi, cmap = plt.cm.Spectral_r, label=r'$w_{f*} = \frac{g \cdot (z_{CL} - z_s) \cdot \Phi \cdot \epsilon }{(\theta_{CL} - \theta_{BL})}}$')
# plt.plot(wStar_zCL, fitZcl[1] + fitZcl[0]*wStar_zCL,c='grey')
# plt.colorbar().set_label('$\Phi$ [Km$^2$/s]')
# plt.legend(fontsize=12)
# ax2.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
# plt.tight_layout()
# plt.savefig(plume.figdir + 'injectionModel/thetaBLinjection.pdf')
# plt.show()
# # plt.close()
#
# # plt.figure()
#


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


#================dimensionless fit=========================
Gamma = np.empty(len(RunList))* np.nan
thetaE =  np.empty(len(RunList))* np.nan
zE = np.empty(len(RunList))* np.nan
for nCase, Case in enumerate(RunList):
    ziidx = np.nanargmin(abs(interpZ - zi[nCase]))
    BLidx = int(ziidx*2/3)
    fittopidx = np.nanargmin(abs(interpZ - 2700))
    baselinedTheta = sounding[nCase,ziidx:fittopidx] - sounding[nCase,BLidx]
    GammaFit = linregress(interpZ[ziidx:fittopidx],baselinedTheta)
    Gamma[nCase] = GammaFit[0]
    thetaE[nCase] = sounding[nCase,BLidx]
    zE[nCase] = -GammaFit[1]/GammaFit[0]
zStar = (zCL - zE)/zi
HStar = (3/4.)**(3/2.)*((thetaE/(g*Gamma**3))**(1/4.)) * np.sqrt((3/2.)*Phi/zi**3)


dimlessFit = linregress(HStar,zStar)
plt.figure()
plt.title('DIMENSIONLESS RELATIONSHIP')
plt.scatter(HStar,zStar,c=Phi,cmap=plt.cm.plasma)
ax = plt.gca()
for i, txt in enumerate(RunList):
    ax.annotate(txt, (HStar[i], zStar[i]),fontsize=6)
plt.gca().set(xlabel = r'$\overline{H}$', ylabel=r'$\overline{z}$')
plt.colorbar(label=r'fireline intensity [Km$^2$/s]')
plt.savefig(plume.figdir + 'injectionModel/DimensionlessGroups.pdf')
plt.show()


#===========iterative solution===============



zCLerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
zCLerrorBiased = np.empty((runCnt)) * np.nan          #parameterization error [m]
zCLcalc = np.empty((runCnt)) * np.nan          #parameterization error [m]

from scipy.optimize import root
for nCase,Case in enumerate(RunList):
    BLfrac = 0.666
    BLidx = np.nanargmin(abs(interpZ - (BLfrac)*zi[nCase]))

    #
    toSolveCase = lambda z : z - (wStarFit8[0]* (BLfrac*zi[nCase]) + wStarFit8[1]) - \
                    wStarFit8[0] * 3/(4*np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaZI[nCase])/(thetaZI[nCase] * (z-zi[nCase]*(2/3.)))))  * \
                    (g*Phi[nCase]*(z-zi[nCase]*(2/3.))*(3/2.)/(thetaZI[nCase] * zi[nCase]))**(1/3.)
    #

    #NOT bias-corrected
    toSolveCaseBiased = lambda z : z - (BLfrac*zi[nCase])  - \
                    3./(4* np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaZI[nCase])/(thetaZI[nCase] * (z-zi[nCase]*(2/3.)))))  * \
                    (g*Phi[nCase]*(z-zi[nCase]*(2/3.))*(3/2.)/(thetaZI[nCase] * zi[nCase]))**(1/3.)

    z_initial_guess = zi[nCase]                   #make initial guess BL height
    z_solution = fsolve(toSolveCase, z_initial_guess,factor=0.1)             #solve
    z_solutionBiased = fsolve(toSolveCaseBiased, z_initial_guess,factor=0.1)             #solve

    zCLcalc[nCase] = (3/4)**6 * (thetaZI[nCase]/g) * (((3/2.)* Phi[nCase]/zi[nCase])**2) * (thetaCL[nCase]-thetaZI[nCase])**(-3) + (2/3.)*zi[nCase]
    # z_solution = root(toSolveCase,z_initial_guess,method='hybr')
    # print(z_solution.x)
    zCLerror[nCase] = z_solution - zCL[nCase]                                #store the solution
    zCLerrorBiased[nCase] = z_solutionBiased - zCL[nCase]                                #store the solution

#

# plt.scatter(zCL,zCLcalc)


plt.figure(figsize=(9,6))

plt.suptitle('ITERATIVE SOLUTION')
gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])


ax0 = plt.subplot(gs[0])
plt.title('Error as f($z_{CL})$: RAW')
plt.scatter(zCL,zCLerrorBiased)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax0.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax1 = plt.subplot(gs[1])
plt.title('Error Statistics')
plt.boxplot(zCLerrorBiased)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])


ax2 = plt.subplot(gs[2])
plt.title('Error as f($z_{CL})$: BIAS CORRECTED')
plt.scatter(zCL,zCLerror)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax2.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax3 = plt.subplot(gs[3])
plt.title('Error Statistics')
plt.boxplot(zCLerror)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax3.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
plt.show()

plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'injectionModel/IterativeSolution.pdf')
plt.show()
plt.close()


# #========sensitivity test for zs values==========================

plt.figure()
ax1 = plt.gca()
ax2 = plt.twinx()
for BLfrac in np.arange(0.5, 0.78, 0.02):
    omega = np.empty((len(RunList)))
    thetas = np.empty((len(RunList)))

    for nCase in range(len(RunList)):
        sidx = np.nanargmin(abs(interpZ - BLfrac*zi[nCase]))
        zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
        omega[nCase] = np.trapz(gradT0interp[nCase,sidx:zCLidx],dx=zstep)
        thetas[nCase] = sounding[nCase,sidx]
    tau = 1/ np.sqrt(g*omega/(thetas * (zCL-zi*BLfrac)))
    wstar =  (3/4) *tau*(g*Phi*(zCL-zi*BLfrac)*(3/2.)/(thetas*zi))**(1/3.)       #attempt at zi_dependent threshold
    wstarfit = linregress(wstar+zi*BLfrac,zCL)
    print(wstarfit)
#
# ziCutoff = 900
# for Zr in range(0,ziCutoff,10):
#     wStartest = (g*Phi[zi>ziCutoff]* (zi[zi>ziCutoff]-Zr)/(Omega[zi>ziCutoff]))**(1/3.)
#     c1, Zr_out, r,p,std = linregress(wStartest,zCL[zi>ziCutoff])
#     print(r)
#     ax1.scatter(Zr,Zr_out,c='C1')
#     ax2.scatter(Zr,r,c='C2')
# l1 = ax1.scatter(Zr,Zr_out,c='C1', label='Zr from fit')
# l2 = ax2.scatter(Zr,r,c='C2', label='R of the fit')
# ax1.set(xlabel='Zr_in [m]',aspect='equal',xlim = [100,ziCutoff],ylim=[100,ziCutoff])
# ax1.set_ylabel('Zr_out [m]', color='C1')
# ax2.set_ylabel('R value',color='C2')
# plt.title('Zi CUTOFF: %s' %ziCutoff)
# plt.legend(handles=[l1,l2])
# plt.show()



#===========gamma solution===============
Gammaerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
GammaerrorBiased = np.empty((runCnt)) * np.nan          #parameterization error [m]
for nCase,Case in enumerate(RunList):
    Gamma_solution = wStarFit8[0]*((3./4)**(3/2.) *(thetaE[nCase]/g)**(1/4.)) * (((3/2.)*Phi[nCase]/zi[nCase])**(0.5)) * (1/Gamma[nCase])**(3/4.) + wStarFit8[0]*zE[nCase] + wStarFit8[1]
    Gammaerror[nCase] = Gamma_solution - zCL[nCase]                                #store the solution
    Gamma_solutionBiased = ((3./4)**(3/2.) *(thetaE[nCase]/g)**(1/4.)) * (((3/2.)*Phi[nCase]/zi[nCase])**(0.5)) * (1/Gamma[nCase])**(3/4.) +zE[nCase]
    GammaerrorBiased[nCase] = Gamma_solutionBiased - zCL[nCase]                                #store the solution


plt.figure(figsize=(9,6))
plt.suptitle('EXPLICIT SOLUTION')
gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
ax0 = plt.subplot(gs[0])
plt.title(r'Error as f($z_{CL})$: RAW')
plt.scatter(zCL,GammaerrorBiased)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax0.set(xlabel=r'$z_{CL}$ [m] ', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax1 = plt.subplot(gs[1])
plt.title('Error Statistics')
plt.boxplot(GammaerrorBiased)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])

ax2 = plt.subplot(gs[2])
plt.title(r'Error as f($z_{CL})$: RAW')
plt.scatter(zCL,Gammaerror)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax2.set(xlabel=r'$z_{CL}$ [m] ', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax3 = plt.subplot(gs[3])
plt.title('Error Statistics')
plt.boxplot(Gammaerror)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax3.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])

plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'injectionModel/ExplicitSolution.pdf')

plt.show()
# plt.close()
