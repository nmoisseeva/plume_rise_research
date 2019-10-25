# nmoisseeva@eoas.ubc.ca
# October 2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

Case = 'W5F7R4'
#=================end of input===============


#----------get preignition temperature profile-----------------
cspath = plume.wrfdir + 'interp/wrfcs_' + Case + '.npy'
print('Opening data: %s' %cspath)
csdict = np.load(cspath, allow_pickle=True).item()

#save initial temperature prfile
profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
profileT = np.mean(csdict['temp'][0,:,:],1)
np.save(profpath,profileT)

#----------check for interpolated data----------------------------
avepath = plume.wrfdir + 'interp/wrfave_' + Case + '.npy'

if os.path.isfile(avepath):
    print('Averaged data found at: %s' %avepath)
    avedict = np.load(avepath,allow_pickle=True).item()   # load here the above pickle
else:
    sys.exit('ERROR: no averaged data found - run prep_plumes.py via submit_interp.sh first on Cedar!')

#mask plume as being at last 30ppm---------------------------------
pm = ma.masked_where(avedict['pm25'] <= 30, avedict['pm25'] )
w = ma.masked_where(avedict['pm25'] <= 30, avedict['w'] )

PMmaxVidx = pm.argmax(0)
xmax,ymax = np.nanargmax(PMmaxVidx), np.nanmax(PMmaxVidx)
centerline = ma.masked_where(plume.lvl[PMmaxVidx] == 0, plume.lvl[PMmaxVidx])

T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

#charatertistics of plume temperature anomalyies---------------------
diffT = ma.masked_where(avedict['pm25'] <= 30, (avedict['temp'].T-T0).T)           #get temperature anomaly
Tctr = np.array([diffT[i,ni] for ni, i in enumerate(PMmaxVidx)])
Wctr = np.array([w[i,ni] for ni, i in enumerate(PMmaxVidx)])
Uctr = np.array([avedict['u'][i,ni] for ni, i in enumerate(PMmaxVidx)])

# #calculate total flux from ground and fire
ignited = np.array([i for i in avedict['ghfx'] if i > 2])

r = len(ignited) * plume.dx

#charageteristic vertical velocity--------------------------------------------
g = 9.81
Fflux = np.mean(ignited) * 1000 / ( 1.2 * 1005)         #get heat flux
# Wf = (g*BLdict['zi'][nCase]*Fflux/BLdict['Ti'][nCase])**(1./3)                     #characteristic velocity based on Deodorff

z1=np.sqrt( (Fflux * g * r**3) / (T0[0] * U0[3]**3) )

# A[nCase] = g* parcelHeat[nCase]/ (zi*Ti[nCase])           #some sort of non-dimensional variable
#===========================plotting===========================
