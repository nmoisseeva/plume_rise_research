#nmoisseeva@eoas.ubc.class
#script to make conserved variable subplots
#May 2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import wrf
from scipy.signal import welch
import imp
from random import randint


#====================INPUT===================
#import all common project variables
import plume
imp.reload(plume) 	#force load each time

Cp= 1005
points = 5000
#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]         #load a list of cases

runCnt = len(RunList)                           #count number of cases

Sd = np.empty((runCnt,points)) * np.nan 				#dry static energy
pm = np.empty((runCnt,points)) * np.nan 				#pm mixing ratio
height = np.empty((runCnt,points)) * np.nan
zi = np.empty((runCnt)) * np.nan

plume_info = np.load(plume.figdir + 'NameZiZcl.npy')
#======================repeat main analysis for all runs first===================
#loop through all LES cases

for nCase,Case in enumerate(RunList):

	#exclude outlier runs that are undealt with
	if Case in plume.exclude_runs:
		continue
	print('Examining case: %s ' %Case)

	endpath = plume.wrfdir + 'interp/end/wrfend_' + Case + '.npy'
	if not os.path.isfile(endpath):
		continue

	print('Opening data: %s' %endpath)
	enddict = np.load(endpath, allow_pickle=True).item()

	#create an interpolated profile of temperature
	if Case[-1:]=='T' or Case[-1:]=='E':
		levels = plume.lvltall
	else:
		levels=plume.lvl

	T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
	zi[nCase] = plume.get_zi(T0)
	dimX = np.shape(enddict['PM25'])[2]

	#RANDOM SAMPLING---------------------------------
	maxPM = np.max(enddict['PM25'],(0,1,2))
	plume_mask = np.where(enddict['PM25'][:,:,:int(dimX/3)] > maxPM*0.01)
	for nTrial in range(points):
		randI = randint(0,len(plume_mask[0])-1)
		nZ,nY,nX = plume_mask[0][randI],plume_mask[1][randI],plume_mask[2][randI]

		Sd[nCase,nTrial] = Cp * (enddict['T'][nZ,nX,nY] + 300)
		pm[nCase,nTrial] = enddict['PM25'][nZ,nX,nY]
		height[nCase,nTrial] = levels[nZ]/zi[nCase]

	plt.figure()
	plt.title('CONSERVED VARIABLE PLOT: %s' %Case)
	plt.scatter(Sd[nCase,:]/1000,pm[nCase,:],c=height[nCase,:], cmap = plt.cm.coolwarm,vmin=0.25,vmax=1.75,s=4)
	plt.gca().set(xlabel='dry static energy ($\Theta \cdot C_p$) [kJ/kg]', ylabel='PM mixing ratio [ug/kg]')
	plt.colorbar(label=r'$z/z_i$')
	# plt.show()
	plt.savefig(plume.figdir + 'mixing/mixing_%s.pdf' %Case)
	plt.close()

	#CENTERLINE SAMPLING---------------------------------
