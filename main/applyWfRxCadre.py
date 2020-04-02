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


#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	             #force load each time
g = 9.81                         #gravity
mf, bf =  6.812, 412.659         #slope and intercept from LES regression model

#=================end of input===============

print('.....extracting NetCDF data from %s ' %plume.rxcadredata)
ncdata = netcdf.netcdf_file(plume.rxcadredata, mode ='r')

#get height data
zstag = (ncdata.variables['PHB'][0,:,:,:] + ncdata.variables['PH'][0,:,:,:])/ g
zdestag = wrf.destagger(zstag,0)
z0 = np.mean(zdestag,(1,2))

#load 10am profile
# sounding_path =
# print('.....extracting vertical profile data from %s' %sounding_path)

#get heat flux and extract only the necessasry portion of the domain
temperature = ncdata.variables['T'][0,:,:,:] + 300.
ghfx = ncdata.variables['GRNHFX'][:,20:70,170:230] 	         #extract fire heat flux
# xtime = ncdata.variables['XTIME'][:] * 60 			#get time in seconds
rotated_fire = rotate(ghfx[:,:,:],125, axes=(1, 2), reshape=False, mode='constant') #near surface wind is at 125deg (atm)


#do fire averaging: use ~6min as interval: estimated as time during which fire remains withing 40m grid given ROS=0.1m/s
masked_flux = ma.masked_less_equal(rotated_fire[70:80,:,:], 500)
ave_masked_flux = np.nanmean(masked_flux,0)
meanFire = np.mean(ave_masked_flux,0)
ignited = np.array([i for i in meanFire if i > 500])       #differs from evalution in RxCADRE where cutoff is 5000 - need to sort this out!!!
Phi = np.trapz(ignited, dx = 40) / ( 1.2 * 1005)     #hardcoded dx=40, to match RxCADRE run

plt.figure()
plt.subplot(121)
plt.imshow(ave_masked_flux,vmin=0, vmax=15000, cmap=plt.cm.gist_heat_r)
plt.subplot(122)
plt.plot(ignited)
plt.show()

#find boundary layer VelScale_InjHeightdef get_zi(T0):
T0 = np.mean(temperature,(1,2))
dT = T0[1:]-T0[0:-1]
dz = z0[1:] - z0[0:-1]
dT_dz = dT/dz
gradT = dT_dz[1:] - dT_dz[0:-1]
si = 6                                                #skipping bottom 100m of the surface inversion
zi_idx = np.argmax(gradT[si:40]) + si                 #vertical level index of BL top
zi = z0[zi_idx]

#interpolate to constant height step
from scipy.interpolate import interp1d
interpz =np.arange(100,2000,10)                     #excluding bottom 100m and interpolating with 10m step
interpf= interp1d(z0[1:], dT)
interpdT = interpf(interpz)


#use numerical solver
# toSolve = lambda z : z - bf - mf * (g*Phi*zi/(np.trapz(dT[si+1:int(z/40.)], dx = 40.)))**(1/3.)
toSolve = lambda z : z - bf - mf * (g*Phi*zi/(np.trapz(interpdT[:int(z/10.)], dx = 10.)))**(1/3.)
z_initial_guess = zi                                        #make initial guess BL height
z_solution = fsolve(toSolve, z_initial_guess)
print(z_solution)
