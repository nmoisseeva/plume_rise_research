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
mf, bf =  6.812, 412.659         #slope and intercept from LES regression model
zstep = 10.                      #vertical step to interpolate to
intFlux = [367, 254, 146, 108, 270, 90 ]    #RxCADRE HIP1 integral flux for HIP1 sensors [kW s /m2]
rawROS = 0.23                    #RxCADRE ROS value for HIP1 (Butler 2016) [m/s]
#=================end of input===============


#get "truth" solution using CWI LES smoke for rxcadre#load and rotate model smoke to mean wind
print('.....extracting interpolated smoke data from %s' %plume.rxinterpCO2)
tracerinterp = np.load(plume.rxinterpCO2)

endingCO2 = np.mean(tracerinterp,0)
singleAngleCO2 = rotate(endingCO2,39, axes=(1, 2), reshape=True, mode='constant')
cwiCO2 = np.mean(singleAngleCO2, 2)
dimZ, dimY, dimX = np.shape(singleAngleCO2)
xx,yy= np.meshgrid(np.arange(0,dimY*40, 40),plume.rxlvl)
plt.title('CWI CO2 AT SIMULATION END')
plt.contourf(xx,yy,cwiCO2,levels = 100, vmax = np.max(cwiCO2)/3, cmap=plt.cm.cubehelix_r)
plt.show()
plt.close()

#locate centerline
ctrZidx = cwiCO2.argmax(0)
ctrXidx = cwiCO2.argmax(1)
co2Ctr = np.array([cwiCO2[ctrZidx[nX],nX] for nX in range(dimY)])
xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
centerline = ma.masked_where(plume.rxlvl[ctrZidx] == 0, plume.rxlvl[ctrZidx])
smoothCenterline = savgol_filter(centerline, 51, 3) # window size 51, polynomial order 3
stableCO2 = cwiCO2[:,120:200]                       #hardcoded based on smoothCenterline
stableProfile = np.mean(stableCO2,1)
rxzCL = np.mean(smoothCenterline[120:200])
rxzCLidx = int(np.mean(ctrZidx[120:200]))

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
ignited = np.array([i for i in meanFire if i > 500])        #differs from evalution in RxCADRE where cutoff is 5000 based on Butler 2016
Phi = np.trapz(ignited, dx = 40.) / ( 1.2 * 1005)           #hardcoded dx=40, to match RxCADRE run

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
# plt.close()

#find boundary layer height
T0 = np.mean(temperature,(1,2))
dTLES = T0[1:]-T0[0:-1]
dz = z0[1:] - z0[0:-1]
dT_dz = dTLES/dz
gradT = dT_dz[1:] - dT_dz[0:-1]
si = 5                                                #skipping bottom 100m of the surface inversion (HARDCODED!)
zi_idx = np.argmax(gradT[si:40]) + si                 #vertical level index of BL top
zi = z0[zi_idx]

#interpolate to constant height step
interpz =np.arange(0,2000,zstep)
interpfLES= interp1d(z0stag[1:-1], dTLES,fill_value='extrapolate')
interpdT_LES = interpfLES(interpz)

plt.figure()
plt.plot(dTLES,z0stag[1:-1])
plt.plot(interpdT_LES,interpz)
plt.show()
plt.close()

#use numerical solver
toSolveLES = lambda z : z - bf - mf * (g*Phi*zi/(np.trapz(interpdT_LES[si:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule
z_initial_guess = zi                                        #make initial guess BL height
LESsolution = fsolve(toSolveLES, z_initial_guess)
print('\033[93m' + 'Smoke injection based on modelled RxCADRE heat and sounding: %.2f' %LESsolution+ '\033[0m')

#---------------------------------SOLUTION USING RAW DATA-----------------------
#import input sounding
#load 10am profile
print('.....extracting vertical profile data from %s' %plume.rxsounding)
sounding = np.genfromtxt(plume.rxsounding, skip_header=1, usecols = [0,1,2,3,4])

#interpolate temperature profile to a constant height step for integration
interpfRAW= interp1d(sounding[:,0], sounding[:,1],fill_value='extrapolate')
interpTRAW = interpfRAW(interpz)
plt.figure()
plt.plot(sounding[:,1],sounding[:,0])
plt.plot(interpTRAW,interpz)
plt.show()
plt.close()

si= int(100/zstep)
dTRAW = interpTRAW[1:] - interpTRAW[:-1]

#get estimates of parameters based on raw data (Phi and zi)
meanQt = np.mean(intFlux,0)
Phi = (meanQt / rawROS) * 4
ziRAW = 1050                       #HARDCODED USING OBSERVATION!!!

#use numerical solver on raw profile----------
toSolveRAW = lambda z : z - bf - mf * (g*Phi*ziRAW/(np.trapz(dTRAW[si:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule
z_initial_guess = ziRAW                                        #make initial guess BL height
RAWsolution = fsolve(toSolveRAW, z_initial_guess)
print('\033[93m' + 'Smoke injection based on raw profile and raw heat: %.f' %RAWsolution + '\033[0m')

#------------------PLOTTING COMPARISON---------------------

plt.figure()
plt.title('EVALUATION OF INJECTION MODEL USING RxCADRE')
plt.plot(stableProfile, plume.rxlvl,label='vertical $CO_2$ profile')
ax = plt.gca()
ax.fill_betweenx(plume.rxlvl, co2Q1, co2Q3, alpha=0.35,label='IQR')
plt.hlines(rxzCL,xmin=0, xmax = 5000,linestyle=':', label='"true" injection height $z_{CL}$ = %.f m' %rxzCL)
plt.hlines(RAWsolution,xmin=0, xmax = 5000,linestyle='--', colors='C3', label='parameterized injection $z_{CL}$ = %.f m' %RAWsolution)
ax.set(xlabel='$CO_2$ concentration [mu/kg]',ylabel='height [m]')
plt.legend(loc='lower right')
plt.savefig(plume.figdir + 'injectionModel/EvaluationRxCADRE.pdf')
plt.show()
