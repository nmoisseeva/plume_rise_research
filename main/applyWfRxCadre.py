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
mf, bf =  2.158, 438.878         #slope and intercept from LES regression model
zstep = 20.                      #vertical step to interpolate to
zs = 200.                       #surface layer height
intFlux = [367, 254, 146, 108, 270, 90 ]    #RxCADRE HIP1 integral flux for HIP1 sensors [kW s /m2]
rawROS = 0.23                    #RxCADRE ROS value for HIP1 (Butler 2016) [m/s]
#=================end of input===============
interpz =np.arange(0,2000,zstep)
si = int(zs/zstep)                                          #skipping surface inversion



#get "truth" solution using CWI LES smoke for rxcadre#load and rotate model smoke to mean wind
print('.....extracting interpolated smoke data from %s' %plume.rxinterpCO2)
tracerinterp = np.load(plume.rxinterpCO2)

endingCO2 = np.mean(tracerinterp,0)                     #average the CO2 over the 10min extracted for modelling
singleAngleCO2 = rotate(endingCO2,39, axes=(1, 2), reshape=True, mode='constant')
cwiCO2 = np.sum(singleAngleCO2, 2) / 1000               #CONVERTED TO PPM!!!! (not ppmv)
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
# rxzCLidx = int(np.mean(ctrZidx[120:200]))
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
ignited = np.array([i for i in meanFire if i > 500])        #differs from evalution in RxCADRE where cutoff is 5000 based on Butler 2016
PhiLES = np.trapz(ignited, dx = 40.) / ( 1.2 * 1005)           #hardcoded dx=40, to match RxCADRE run

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
interpT0LES = interpfLES(interpz)

gradTLES = (interpT0LES[1:] - interpT0LES[:-1])/zstep
gradT2 = gradTLES[1:] - gradTLES[:-1]
zi_idx = np.argmax(gradT2[si:]) + si
zi = interpz[zi_idx]

plt.figure()
plt.plot(T0,z0)
plt.plot(interpT0LES,interpz)
plt.show()
plt.close()

#use numerical solver
toSolveLES = lambda z : z - bf - mf * (g*PhiLES*(zi-zs)/(np.trapz(gradTLES[si:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule
# toSolveLES = lambda z : z - bf - mf * (g*PhiLES*(zi-zs)/(interpT0LES[int(z/zstep)] - interpT0LES[si]))**(1/3.)           #using delta

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

gradTRAW = (interpTRAW[1:] - interpTRAW[:-1])/zstep

#get estimates of parameters based on raw data (Phi and zi)
meanQt = np.mean(intFlux,0)
PhiRAW = (meanQt / rawROS) * 4
ziRAW = 1050                       #HARDCODED USING OBSERVATION!!!

#use numerical solver on raw profile----------
toSolveRAW = lambda z : z - bf - mf * (g*PhiRAW*(ziRAW-zs)/(np.trapz(gradTRAW[si+1:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule
# toSolveRAW = lambda z : z - bf - mf * (g*PhiRAW*(ziRAW-zs)/(interpTRAW[int(z/zstep)]-interpTRAW[si]))**(1/3.)           #using delta

z_initial_guess = ziRAW                                        #make initial guess BL height
RAWsolution = fsolve(toSolveRAW, z_initial_guess)
print('\033[93m' + 'Smoke injection based on raw profile and raw heat: %.f' %RAWsolution + '\033[0m')


#-----------------solve with RAW heat and LES T0 and CWI Smoke

#use numerical solver
toSolveLESRAW = lambda z : z - bf - mf * (g*PhiRAW*(zi-zs)/(np.trapz(gradTLES[si:int(z/zstep)], dx = zstep)))**(1/3.)           #using trapezoidal rule

z_initial_guess = zi                                        #make initial guess BL height
LESRAWsolution = fsolve(toSolveLESRAW, z_initial_guess)
print('\033[93m' + 'Smoke injection based on modelled RxCADRE heat and sounding: %.2f' %LESRAWsolution+ '\033[0m')

#------------------PLOTTING COMPARISON---------------------

plt.figure()
plt.title('EVALUATION OF INJECTION MODEL USING RxCADRE')
plt.plot(stableProfile, plume.rxlvl,label='vertical $CO_2$ profile')
ax = plt.gca()
ax.fill_betweenx(plume.rxlvl, co2Q1, co2Q3, alpha=0.35,label='IQR')
plt.hlines(rxzCL,xmin=0, xmax = 3000, label='"true" injection height $z_{CL}$ = %.f m' %rxzCL)
plt.hlines(RAWsolution,xmin=0, xmax = 3000,linestyle='--', colors='C1', label='RAW (heat and sounding) $z_{CL}$ = %.f m' %RAWsolution)
plt.hlines(LESsolution,xmin=0, xmax = 3000,linestyle=':', colors='C2', label='LES (heat and sounding) $z_{CL}$ = %.f m' %LESsolution)
plt.hlines(LESRAWsolution,xmin=0, xmax = 3000,linestyle='--', colors='C3', label='LES (sounding) and RAW (heat) $z_{CL}$ = %.f m' %LESRAWsolution)
ax.set(xlabel='$CO_2$ concentration [ppm]',ylabel='height [m]')
plt.legend(loc='lower right', fontsize=9.5)
plt.savefig(plume.figdir + 'injectionModel/EvaluationRxCADRE.pdf')
plt.show()
