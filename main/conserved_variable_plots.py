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
points = 1000
#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]         #load a list of cases
runCnt = len(RunList)                           #count number of cases

Sd = np.empty((runCnt,points)) * np.nan 				#dry static energy
pm = np.empty((runCnt,points)) * np.nan 				#pm mixing ratio
height = np.empty((runCnt,points)) * np.nan

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

	#get dimension of domain and pick a random sample
	dimY, dimX = np.shape(enddict['GRNHFX'])
	for nTrial in range(points):
		randY, randX = randint(int(dimY/4), 3*int(dimY/4)), randint(0, int(dimX/4))
		randZ = randint(np.argmin(abs(plume.lvl-200)), np.argmin(abs(plume.lvl-1000)))

		Sd[nCase,nTrial] = Cp * (enddict['T'][randZ,randY,randX] + 300)
		pm[nCase,nTrial] = enddict['PM25'][randZ,randY,randX]
		height[nCase,nTrial] = plume.lvl[randZ]

	plt.figure()
	plt.title('CONSERVED VARIABLE PLOT: %s' %Case)
	plt.scatter(Sd[nCase,:]/1000,pm[nCase,:],c=height[nCase,:])
	plt.gca().set(xlabel='dry static energy ($\Theta \cdot C_p$) [kJ/kg]', ylabel='PM mixing ratio [ug/kg]')
	plt.colorbar(label='height AGL [m]')
	# plt.show()
	plt.savefig(plume.figdir + 'mixing/mixing_%s.pdf' %Case)
	plt.close()
