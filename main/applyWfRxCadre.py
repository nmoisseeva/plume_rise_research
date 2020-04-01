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


#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

#=================end of input===============

print('.....extracting NetCDF data from %s ' %plume.rxcadredata)
ncdata = netcdf.netcdf_file(plume.rxcadredata, mode ='r')

#get heat flux and extract only the necessasry portion of the domain
ghfx = ncdata.variables['GRNHFX'][:,20:70,170:230] 	         #extract fire heat flux
# xtime = ncdata.variables['XTIME'][:] * 60 			#get time in seconds

rotated_fire = rotate(ghfx[:,:,:],125, axes=(1, 2), reshape=False, mode='constant') and #near surface wind is at 125deg (atm)

#plot snapshot of rotated fire
plt.imshow(rotated_fire[70,:,:], vmin=0, vmax=15000, cmap=plt.cm.gist_heat_r)
plt.colorbar()
plt.show()

#do fire averaging: use ~6min as interval: estimated as time during which fire remains withing 40m grid given ROS=0.1m/s
masked_flux = ma.masked_less_equal(rotated_fire[70:100,:,:], 0)
ave_masked_flux = np.nanmean(masked_flux,0)
# plt.imshow(ave_masked_flux,vmin=0, vmax=15000, cmap=plt.cm.gist_heat_r)
cs_flux = np.mean(ave_masked_flux,0)




print('..... creating an igntion mask on atm grid (may take several minutes)')
ign_mask_atm = np.empty_like(ghfx) * np.nan
for nt in range(len(xtime)):
	print(nt)
	current_ign = ghfx[nt,:,:]
	temp_mask = np.empty_like(current_ign) * np.nan
	temp_mask[current_ign>5000] = 1 	#residence time defined as at least 5kW/m2 as per Butler2013
	ign_mask_atm[nt,:,:] = temp_mask

    print('WARNING: Slow routine: rotating the array to be alighned with mean wind')


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

# #rotating the array
# cw_sum = []
# for nTime in range(nT):
#     print(nTime)
#     CO2rot = rotate(tracerinterp[nTime,:,:,:], 39, axes=(1, 2), reshape=True, mode='constant')
#     CO2tot = np.sum(CO2rot,2)
#     cw_sum.append(CO2tot)
#
# cw_sum = np.array(cw_sum)
# fig = plt.figure()
# ax = plt.gca()
