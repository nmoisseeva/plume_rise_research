# April 2020
#nmoisseeva@eoas.ubc.ca
#This code applies plume injection model to RxCADRE data


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.optimize import fsolve
from matplotlib import gridspec
from scipy.ndimage.interpolation import rotate
import wrf
import matplotlib.dates as mdates
import datetime as dt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	             #force load each time
g = 9.81                         #gravity
# mf, bf =  2.158, 438.878         #slope and intercept from LES regression model
mf, bf = 0.9243, 115.059
C = 1.0087
zstep = 20.                      #vertical step to interpolate to
# zs = 200.                       #surface layer height
BLfrac = 0.75           #fraction of BL height to set zs at
# intFlux = [367, 254, 146, 348, 270, 90]    #RxCADRE HIP1 integral flux for HIP1 sensors [kW s /m2] FBP ID [22,2,20,3,14,19]
intFlux = [367,254,146,348,    357,258,256,311,126,    279,225,295,284,472]    #RxCADRE HIP1 integral flux for HIP1 sensors [kW s /m2] FBP ID [22,2,20,3,   5,13,15,7,4,   12,10,16,11,8]
rawROS = np.mean([0.32,0.24,0.24,0.31,0.25,0.44,0.22,0.61,     0.47,0.38,0.47,0.47,0.44,0.45,0.48,    0.25,0.46,0.39,0.46,0.16,0.51,0.3,0.38,0.27,0.43])

# rawROS = np.mean([0.23,0.40])                    #RxCADRE ROS value for HIP1 (Butler 2016) [m/s]
rawROS = np.mean([0.23,0.4,0.44,0.36,0.23,0.42])
#=================end of input===============
interpz =np.arange(0,2000,zstep)
# si = int(zs/zstep)                                          #skipping surface inversion

#get "truth" solution using CWI LES smoke for rxcadre#load and rotate model smoke to mean wind
print('.....extracting interpolated smoke data from %s' %plume.rxinterpCO2)
tracerinterp = np.load(plume.rxinterpCO2)

# endingCO2 = np.mean(tracerinterp,0)                     #average the CO2 over the 10min extracted for modelling
endingCO2 = tracerinterp[0,:,:,:]                   #average the CO2 over the 10min extracted for modelling

singleAngleCO2 = rotate(endingCO2,39, axes=(1, 2), reshape=True, mode='constant')
cwiCO2 = np.sum(singleAngleCO2, 2) / 1000               #CONVERTED TO PPM!!!! (not ppmv)
dimZ, dimY, dimX = np.shape(singleAngleCO2)
xx,yy= np.meshgrid(np.arange(0,dimY*40, 40),plume.rxlvl)


cleanCO2 = ma.masked_where(cwiCO2 <= 1000, cwiCO2 ) #mask all non-plume cells
#locate centerline
ctrZidx = cleanCO2.argmax(0)
ctrXidx = cleanCO2.argmax(1)
co2Ctr = np.array([cwiCO2[ctrZidx[nX],nX] for nX in range(dimY)])
xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
centerline = ma.masked_where(plume.rxlvl[ctrZidx] == 0, plume.rxlvl[ctrZidx])

# filter_window = max(int(plume.read_tag('W',[Case])*10+1),51)
smoothCenterline = savgol_filter(centerline, 51, 3)

#this portion doesn't work - needs debugging------------------
#calculate concentration changes along the centerline
dCO2dX = co2Ctr[1:]-co2Ctr[0:-1]
smoothCO2 = savgol_filter(dCO2dX, 51, 3) # window size 101, polynomial order 3
stableCO2mask = [True if abs(smoothCO2[nX])< np.nanmax(smoothCO2)*0.1 and \
                        abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                        nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                        nX > np.nanargmax(smoothCO2) and\
                        nX > np.nanargmax(centerline) +10 and
                        nX > np.nanargmax(smoothCenterline)+10 else \
                        False for nX in range(dimY-1) ]
stableCO2 = cleanCO2[:,1:][:,stableCO2mask]
stableProfile = np.mean(stableCO2,1)
#----------------- needs debugging------------------

stableCO2 = cwiCO2[:,130:160]                       #hardcoded based on smoothCenterline
stableProfile = np.mean(stableCO2,1)



#plot as sanity check
plt.figure(figsize=(15,4))
plt.title('CWI CO2 AT SIMULATION END')
plt.contourf(xx,yy,cwiCO2,levels = 100, vmax = np.max(cwiCO2)/3, cmap=plt.cm.cubehelix_r)
plt.plot(xx[0],centerline,':')
plt.plot(xx[0], smoothCenterline)
ax = plt.gca()
ax.set(aspect='equal')
# ax.fill_between(xx[0][1:], 0, 1, where=stableCO2mask, color='grey', alpha=0.4, transform=ax.get_xaxis_transform(), label='averaging window')
plt.tight_layout()
plt.show()
plt.close()


rxzCL = np.mean(smoothCenterline[130:160])
rxzCLidx = np.argmin(abs(interpz - rxzCL))

#get stats and plot profile
co2Q1 = np.percentile(stableCO2,25,axis = 1)
co2Q3 = np.percentile(stableCO2,75,axis = 1)

print('\033[93m' + '"True" injection based on LES CWI smoke: %.2f' %rxzCL+ '\033[0m')

#-------------------------------------SOLUTION USING LES SIMULATION-----------------
print('.....extracting NetCDF data from %s ' %plume.rxcadredata)
ncdata = netcdf.netcdf_file(plume.rxcadredata, mode ='r')

#get height data
zstag = (ncdata.variables['PHB'][0,:,:,:] + ncdata.variables['PH'][0,:,:,:])/ g
zdestag = wrf.destagger(zstag,0)
z0 = np.mean(zdestag,(1,2))
z0stag = np.mean(zstag,(1,2))

#get heat flux and extract only the necessasry portion of the domain
temperature = ncdata.variables['T'][0,:,:,:] + 300.
ghfx = ncdata.variables['GRNHFX'][:,20:70,170:230] 	         #extract fire heat flux
rotated_fire = rotate(ghfx[:,:,:],125, axes=(1, 2), reshape=False, mode='constant') #near surface wind is at 125deg (atm)

#do fire averaging: use ~6min as interval: estimated as time during which fire remains withing 40m grid given ROS=0.1m/s
masked_flux = ma.masked_less_equal(rotated_fire[60:80,:,:], 500)
ave_masked_flux = np.nanmean(masked_flux,0)
meanFire = np.mean(ave_masked_flux,0)
ignited = np.array([i for i in meanFire if i > 1000])        #differs from evalution in RxCADRE where cutoff is 5000 based on Butler 2016
PhiLES = np.trapz(ignited, dx = 40.) / (1.2 * 1005)           #hardcoded dx=40, to match RxCADRE run

plt.figure(figsize=(12,8))
plt.subplot(121)
plt.title('HEAT FLUX INPUT')
plt.imshow(ave_masked_flux,vmin=0, vmax=15000, cmap=plt.cm.gist_heat_r)
plt.subplot(122)
plt.title('HEAT FLUX CROSSSECTION')
plt.plot(ignited)
plt.gca().set(xlabel='grid number', ylabel='mean heat flux [W/m2]')
plt.tight_layout()
plt.show()
plt.close()

#find boundary layer height
T0 = np.mean(temperature,(1,2))
interpfLES= interp1d(z0, T0,fill_value='extrapolate')
soundingLES = interpfLES(interpz)

gradTLES = (soundingLES[1:] - soundingLES[:-1])/zstep
gradT2 = gradTLES[1:] - gradTLES[:-1]
ziidx = np.argmax(gradT2[10:]) + 10
zi = interpz[ziidx]
zS = zi*BLfrac
zsidx = np.argmin(abs(interpz - zS))

plt.figure()
plt.plot(T0,z0)
plt.plot(soundingLES,interpz)
plt.show()
plt.close()

thetaS = soundingLES[zsidx]
thetaCL = soundingLES[rxzCLidx]

#use numerical solver
# toSolveLES = lambda z : z - (mf*BLfrac*zi + bf) - \
#                     mf * 3/(4*np.sqrt(g*(soundingLES[int(z/zstep)] - thetaS)/(thetaS * (z-zi*BLfrac))))  * \
#                     (g*PhiLES*(z-zi*BLfrac)*(3/2.)/(thetaS * zi))**(1/3.)

toSolveLES = lambda z : z - bf - mf* (zS + \
                    C/(np.sqrt(g*(soundingLES[int(z/zstep)] - thetaS)/(thetaS * (z-zS))))  * \
                    (g*PhiLES*(z-zS)/(thetaS * zi))**(1/3.))

z_initial_guess = zi                                        #make initial guess BL height
LESsolution = fsolve(toSolveLES, z_initial_guess,factor = 0.1)
print('\033[93m' + 'Smoke injection based on modelled RxCADRE heat and sounding: %.2f' %LESsolution+ '\033[0m')

#-----------------solve with RAW heat and LES T0 and CWI Smoke
#get estimates of parameters based on raw data (Phi and zi)
meanQt = np.mean(intFlux,0) * 1000 / (1005 * 1.2)
PhiRAW = round(meanQt) * round(rawROS,2) * 4

#use numerical solver
# toSolveLESRAW = lambda z : z - bf - mf * (g*PhiRAW*(zi-zs)/(np.trapz(gradTLES[si:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule

toSolveLESRAW = lambda z : z - bf - mf* (zS + \
                    C/(np.sqrt(g*(soundingLES[int(z/zstep)] - thetaS)/(thetaS * (z-zS))))  * \
                    (g*PhiRAW*(z-zS)/(thetaS * zi))**(1/3.))

z_initial_guess = zi                                        #make initial guess BL height
LESRAWsolution = fsolve(toSolveLESRAW, z_initial_guess,factor=0.1)
print('\033[93m' + 'Smoke injection based on modelled RxCADRE heat and sounding: %.2f' %LESRAWsolution+ '\033[0m')

#------------------PLOTTING COMPARISON---------------------

plt.figure()
plt.title('EVALUATION OF INJECTION MODEL USING RxCADRE')
plt.plot(stableProfile, plume.rxlvl,label='vertical $CO_2$ profile')
ax = plt.gca()
ax.fill_betweenx(plume.rxlvl, co2Q1, co2Q3, alpha=0.35,label='IQR')
# plt.hlines(rxzCL,xmin=0, xmax = 3000, label='"true" injection height $z_{CL}$ = %.f m' %rxzCL)
# plt.hlines(LESsolution,xmin=0, xmax = 3000,linestyle=':', colors='C2', label='LES (sounding and heat) $z_{CL}$ = %.f m' %LESsolution)
# plt.hlines(LESRAWsolution,xmin=0, xmax = 3000,linestyle='--', colors='C3', label='LES (sounding) and RAW (heat) $z_{CL}$ = %.f m' %LESRAWsolution)
plt.hlines(rxzCL,xmin=0, xmax = 3000, label='LES $z_{CL}$ = %.f m' %round(rxzCL,-1))
plt.hlines(LESsolution,xmin=0, xmax = 3000,linestyle='solid', colors='C1', label='modelled $z_{CL}$ = %.f m,  using $I_{LES}$' %round(LESsolution[0],-1))
plt.hlines(LESRAWsolution,xmin=0, xmax = 3000,linestyle='dashed', colors='C3', label='modelled $z_{CL}$ = %.f m,  using $I_{obs}$' %round(LESRAWsolution[0],-1))
ax.set(xlabel='$CO_2$ concentration [ppm]',ylabel='height [m]')
plt.legend(loc='lower right', fontsize=9.5)
plt.savefig(plume.figdir + 'injectionModel/EvaluationRxCADRE.pdf')
plt.show()
