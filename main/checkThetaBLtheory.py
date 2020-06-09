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
Omega = np.empty((runCnt)) * np.nan             #cumulative vertical temperature (Km) - the questionable denominator term
sounding = np.empty((runCnt,len(interpZ))) * np.nan         #storage for interpolated soundings
gradT0interp = np.empty((runCnt,len(interpZ)-1)) * np.nan   #storage for temperature gradient
wZi = np.empty((runCnt)) * np.nan
wCum = np.empty((runCnt)) * np.nan
wCum2 = np.empty((runCnt)) * np.nan
thetaFZi = np.empty((runCnt)) * np.nan
thetaZi = np.empty((runCnt)) * np.nan
r = np.empty((runCnt)) * np.nan

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
    r[nCase] = len(ignited) * plume.dx
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
    Omega[nCase] = np.trapz(dTinterp[siBL:zCLidx], dx = zstep)     #calculate the denominator term by integrating temperature change with height, excluding surface layer (Km)

    w = csdict['w'][-1,:,:]
    Wctr = np.array([w[nZ, ctrXidx[nZ]] for nZ in range(dimZ)])    #get concentration along the centerline
    # Wcum[nCase] = np.sum(np.max(w,1)[:ymax+5])
    if Case[-1:]=='T' or Case[-1:]=='E':
        interpWctr = interp1d(plume.lvltall, Wctr,fill_value='extrapolate')
    else:
        interpWctr= interp1d(plume.lvl, Wctr,fill_value='extrapolate')
    Wctrinterp = interpWctr(interpZ)
    wZi[nCase] = np.sum(Wctrinterp[ziidx])
    wCum[nCase] = np.sum(np.max(w,1)[25:zCLidx])
    wCum2[nCase] = np.sum(Wctrinterp[25:zCLidx])
    #do theta testing
    temperature = csdict['temp'][-1,:,:]
    Tctr = np.array([temperature[nZ, ctrXidx[nZ]] for nZ in range(dimZ)])    #get concentration along the centerline

    if Case[-1:]=='T' or Case[-1:]=='E':
        interpTctr = interp1d(plume.lvltall, Tctr,fill_value='extrapolate')
    else:
        interpTctr= interp1d(plume.lvl, Tctr,fill_value='extrapolate')
    Tctrinterp = interpTctr(interpZ)
    thetaFZi[nCase] = Tctrinterp[ziidx]
    thetaZi[nCase] = sounding[nCase,ziidx]
#======================compare model formulations========================
# plt.hist(thetaFZi-thetaZi,bins=20,alpha=0.5)
# plt.hist(Omega,bins=20,alpha=0.5)
# plt.show()
plt.figure()
plt.hist(wCum-wCum2,bins=20)
plt.show()
# plt.close()

#define wf* (as per original 'wrong' formulation)
wStar = (g*Phi*(zi-200)*(3./2)/(Omega))**(1/3.)

wZicalc = (g*Phi*zi*(3./2)/(thetaZi*r))**(1/3.)

plt.figure()
plt.scatter(wZi,wZicalc)
plt.gca().set(aspect='equal')


# plt.figure()
# plt.scatter(wStar,wCum)
# plt.gca().set(aspect='equal')
# plt.show()

# plt.scatter(wStar,Wcum)
# plt.plot(np.arange(1000),np.arange(1000))
# plt.show()
# plt.close()

#do linear regression using all data
slopeALL, interceptALL, r_valueALL, p_valueALL, std_errALL = linregress(wStar[np.isfinite(wStar)],zCL[np.isfinite(wStar)])


wStar_zi = (g*Phi*(zi-200)*(3./2)/(Omega))**(1/3.)
wStar_zCL = (g*Phi*(zCL-200)*(3./2)/(Omega))**(1/3.)
#do linear regression using all data
fitZi = linregress(wStar_zi,zCL)
fitZcl = linregress(wStar_zCL,zCL)


#make scatterplot comparisons
plt.figure(figsize=(12,5))
plt.subplot(121)
ax1 = plt.gca()
plt.title('Wf* using zi: [R=%0.2f]' %fitZi[2])
plt.scatter(wStar_zi, zCL, c=Phi, cmap = plt.cm.Spectral_r, label=r'$w_{f*} = \frac{g \cdot (z_i - z_s) \cdot \Phi \cdot \epsilon }{(\theta_{CL} - \theta_{BL})}}}$')
plt.plot(wStar_zi, fitZi[1] + fitZi[0]*wStar_zi,c='grey')
plt.colorbar().set_label('$\Phi$ [Km$^2$/s]')
plt.legend(fontsize=12)
ax1.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
plt.subplot(122)
ax2 = plt.gca()
plt.title('Wf* using zCL: [R=%0.2f]' %fitZi[2])
plt.scatter(wStar_zCL, zCL, c=Phi, cmap = plt.cm.Spectral_r, label=r'$w_{f*} = \frac{g \cdot (z_{CL} - z_s) \cdot \Phi \cdot \epsilon }{(\theta_{CL} - \theta_{BL})}}$')
plt.plot(wStar_zCL, fitZcl[1] + fitZcl[0]*wStar_zCL,c='grey')
plt.colorbar().set_label('$\Phi$ [Km$^2$/s]')
plt.legend(fontsize=12)
ax2.set(xlabel='$w_{f*}$ [m/s]',ylabel='zCL [m]')
plt.tight_layout()
plt.savefig(plume.figdir + 'injectionModel/thetaBLinjection.pdf')
plt.show()
# plt.close()

# plt.figure()

#======================train and test regression model===================
