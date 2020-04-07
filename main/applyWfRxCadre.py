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

#=================end of input===============
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
ignited = np.array([i for i in meanFire if i > 500])       #differs from evalution in RxCADRE where cutoff is 5000 - need to sort this out!!!
Phi = np.trapz(ignited, dx = 40.) / ( 1.2 * 1005)     #hardcoded dx=40, to match RxCADRE run

plt.figure()
plt.subplot(121)
plt.imshow(ave_masked_flux,vmin=0, vmax=15000, cmap=plt.cm.gist_heat_r)
plt.subplot(122)
plt.plot(ignited)
plt.show()
# plt.close()

#find boundary layer
T0 = np.mean(temperature,(1,2))
dT = T0[1:]-T0[0:-1]
dz = z0[1:] - z0[0:-1]
dT_dz = dT/dz
gradT = dT_dz[1:] - dT_dz[0:-1]
si = 5                                                #skipping bottom 100m of the surface inversion
zi_idx = np.argmax(gradT[si:40]) + si                 #vertical level index of BL top
zi = z0[zi_idx]

#interpolate to constant height step
from scipy.interpolate import interp1d
interpz =np.arange(0,2000,10)
interpf= interp1d(z0stag[1:-1], dT,fill_value='extrapolate')
interpdT = interpf(interpz)
plt.plot(dT,z0stag[1:-1])
plt.plot(interpdT,interpz)
plt.show()
plt.close()

#use numerical solver
toSolve = lambda z : z - bf - mf * (g*Phi*zi/(np.trapz(interpdT[si:int(z/10.)], dx = 10.)))**(1/3.)           #using trapezoidal rule
# toSolve = lambda z : z - bf - mf * (g*Phi*zi/(sum(dT[si:(np.abs(z0-z)).argmin()]*dz[si:(np.abs(z0-z)).argmin()])))**(1/3.)             #on model levels using sum

z_initial_guess = zi                                        #make initial guess BL height
z_solution = fsolve(toSolve, z_initial_guess)
print(z_solution)

#---------------------------------SOLUTION USING RAW DATA-----------------------
basetime = dt.datetime(year=2012,month=11,day=10)
#extract and format dispersion data
print('.....importing dispersion data from %s' %plume.rxdispersion)
disp_dict = {}
disp_array = np.genfromtxt(plume.rxdispersion, skip_header=1, usecols = [1,2,3,4,5,7,8,9], delimiter=',')
disp_dict['time'] = np.array([basetime + dt.timedelta(seconds = int(i)) for i in disp_array[:,0]])
disp_dict['CO'] = disp_array[:,1]
disp_dict['CO2'] = disp_array[:,2]
disp_dict['CH4'] = disp_array[:,3]
disp_dict['H2O'] = disp_array[:,4]
disp_dict['lcn'] = np.array(list(zip(disp_array[:,5],disp_array[:,6],disp_array[:,7]-plume.rxsfchgt)))
disp_dict['meta']= 'time: min since restart run | \
					CO: Mixing ratio of carbon monoxide in units of parts per million by volume (ppmv) in dry air. | \
					CO2: Mixing ratio of carbon dioxide in units of ppmv in dry air. | \
					CH4: Mixing ratio of methane in units of ppmv in dry air. | \
					H2O: Mixing ratio of water vapor in percent by volume. | \
					lcn: (lat, lon, elevation) - coords in WGS84, elevation AGL'

plt.figure()
plt.title('AIRPLANE HEIGHT')
plt.plot(disp_dict['time'],disp_dict['lcn'][:,2])
ax = plt.gca()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
plt.gcf().autofmt_xdate()
plt.xlabel('time (CST)')
plt.ylabel('airplane height [m]')
plt.ylim([0,2000])
# plt.savefig(rx.fig_dir + 'PlaneHeightProfiles.pdf')
plt.show()
plt.close()



#import input sounding
#load 10am profile
print('.....extracting vertical profile data from %s' %plume.rxsounding)
sounding = np.genfromtxt(plume.rxsounding, skip_header=1, usecols = [0,1,2,3,4])

#rotate model smoke to mean wind----------
endingCO2 = np.mean(ncdata.variables['tr17_1'][250:,:,:,:],0)
singleAngleCO2 = rotate(endingCO2,39, axes=(1, 2), reshape=True, mode='constant')
cwiCO2 = np.mean(singleAngleCO2, 2)
dimZ, dimY, dimX = np.shape(singleAngleCO2)
xx,yy= np.meshgrid(np.arange(0,dimY*40, 40),z0)
plt.title('CWI CO2 AT SIMULATION END')
plt.contourf(xx,yy,cwiCO2,levels = 100, vmax = np.max(cwiCO2)/2, cmap=plt.cm.cubehelix_r)
plt.show()

#locate centerline
ctrZidx = cwiCO2.argmax(0)
ctrXidx = cwiCO2.argmax(1)
co2Ctr = np.array([cwiCO2[ctrZidx[nX],nX] for nX in range(dimY)])
for nZ in range(dimZ):
    if co2Ctr[ctrXidx[nZ]] < cwiCO2[nZ,ctrXidx[nZ]]:
        co2Ctr[ctrXidx[nZ]] = cwiCO2[nZ,ctrXidx[nZ]]
xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)
centerline = ma.masked_where(z0[ctrZidx] == 0, z0[ctrZidx])
smoothCenterline = savgol_filter(centerline, 51, 3) # window size 51, polynomial order 3

stableCO2 = cwiCO2[:,150:250]                       #hardcoded based on smoothCenterline
stableProfile = np.mean(stableCO2,1)

rxzCL = np.mean(smoothCenterline[150:250])
rxzCLidx = int(np.mean(ctrZidx[150:250]))

co2Q1 = np.percentile(stableCO2,25,axis = 1)
co2Q3 = np.percentile(stableCO2,75,axis = 1)
plt.plot(stableProfile, z0,label='vertical CO2 profile')
plt.hlines(rxzCL,xmin=0, xmax = dimY*40,linestyle=':', label='injection based on rxcadre')
plt.gca().fill_betweenx(z0, co2Q1, co2Q3, alpha=0.35,label='IQR')
plt.legend()
plt.show()


print('Smoke injection based on modelled RxCADRE: %.2f' %rxzCL)
#
#
#
#
# #HARDCODED:
# rotCS_lcn = [30.56653,-86.75507]    #calculated by rotating 30
# ROT_cs_dist, ROT_cs_grid_id = gridTree.query(np.array(rotCS_lcn))
# ROTidxy,ROTidxx = np.unravel_index(ROT_cs_grid_id,np.shape(wrfgeo['XLONG']))
# ROT_profile_mukg = tracerinterp[258,:,ROTidxy,ROTidxx]
# ROT_profile = 1e-3 * ROT_profile_mukg * molar_mass['air'] / molar_mass['CO2']    #convert to ppmv, assuming model is in mug/kg
# #------end hardcoding
#
# cs_lon =  np.mean(disp_dict['lcn'][tidx][ics:fcs,1])
# cs_lat = np.mean(disp_dict['lcn'][tidx][ics:fcs,0])
# tot_column_mukg = np.sum(ncdict['CO2'][258,:,:,:],0)
# smokeim = 1e-3 * tot_column_mukg * molar_mass['air']/molar_mass['CO2']   #convert to ppmv
# im = bm.imshow(smokeim, cmap = plt.cm.bone_r, origin='lower',vmin=0,vmax=801)
# bm.scatter(cs_lon,cs_lat,40,marker='*',color='r',label='average location of corkscrew profile')
# bm.scatter(rotCS_lcn[1],rotCS_lcn[0],20,marker='o',color='k',label='wind-corrected profile location')
# plt.colorbar(im, label='total column CO$_2$ (ppmv)')
# plt.legend()
# plt.tight_layout()
# plt.savefig(rx.fig_dir + 'CSLocation.pdf')
# plt.show()
# plt.close()


#interpolate temperature profile to a constant height step for integration
interpf= interp1d(sounding[:,0], sounding[:,1],fill_value='extrapolate')
interpT = interpf(interpz)
plt.plot(sounding[:,1],sounding[:,0])
plt.plot(interpT,interpz)
plt.show()
plt.close()
si=3
dT = interpT[1:] - interpT[:-1]

# #raw data Phi
intFlux = [367, 254, 146, 108, 270, 90 ] #all values from integrated total flux (in kW s /m2) from HIP1
rawROS = 0.25
meanQt = np.mean(intFlux,0)
Phi = (meanQt / rawROS) * 4

#use numerical solver on raw profile----------
toSolve = lambda z : z - bf - mf * (g*Phi*zi/(np.trapz(dT[si:int(z/10.)], dx = 10.)))**(1/3.)           #using trapezoidal rule
z_initial_guess = zi                                        #make initial guess BL height
z_solution = fsolve(toSolve, z_initial_guess)
print('Smoke injection based on raw profile and modelled fire: %.2f' %z_solution)
