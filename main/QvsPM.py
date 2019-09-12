# nmoisseeva@eoas.ubc.ca
# June 2018


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import sys
import imp
import wrf
#====================INPUT===================
#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

testpath = plume.wrfdir + 'wrfout_W3S200F7R2_short'
#=================end of input===============

testdata = netcdf.netcdf_file(testpath, mode ='r')
ncdict = wrf.extract_vars(testdata, None, ('QVAPOR','tr17_1'))
ncdict['PM25'] = ncdict.pop('tr17_1')

plt.plot(ncdict['PM25'][70,:,75,100],'r--')
ax1 = plt.gca()
ax1.set_ylim([0,max(ncdict['PM25'][70,:,75,100])])
ax2 = plt.gca().twinx()
plt.plot(ncdict['QVAPOR'][70,:,75,100])
ax2.set_ylim([0,max(ncdict['QVAPOR'][70,:,75,100])])
plt.show()
