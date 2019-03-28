#nmoisseeva@eoas.ubc.ca
#March 2019

import numpy as np
import matplotlib.pyplot as plt
import os.path
import wrf
import imp
import plume        #all our data


#====================INPUT===================
#all common variables are stored separately
imp.reload(plume) 	#force load each time
#=================end of input===============

print('ANALYSIS OF PLUME CROSS-SECTIONS')
print('===================================')

for nCase,Case in enumerate(plume.tag):
	print('Examining case: %s ' %Case)

	#----------check for interpolated data----------------------------
	interppath = plume.wrfdir + 'interp/wrfinterp_' + Case + '.npy'
    wrfpath = wrfdir + 'wrfout_'+ tag
    wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')
    ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX'))

	if os.path.isfile(interppath):
		print('Interpolated data found at: %s' %interppath)
		interpdict = np.load(interppath).item()   # load here the above pickle
	else:
		sys.exit('ERROR: no interpolated data found - run prep_plumes.py first!')

    #convert and average data------------------------------------------
    ghfx = ncdict['GRNHFX']/1000. 			#convert to kW
    qvapor = interpdict['QVAPOR']*1000.		#convert to g/kg
    temp = interpdict['T']+300. 			#add perturbation and base temperature
    w = interpdict['W']

    #get dimensions
    dimt, dimy, dimx = np.shape(ghfx)
    xsx = int(round(dimy/2.))

    #create fire cross-section averages
    ghfx_mave = np.mean(ghfx[:,xsx-plume.cs:xsx+plume.cs,:],1)
    w_mave = np.mean(w[:,:,xsx-plume.cs:xsx+plume.cs,:],2)
    temp_mave = np.mean(temp[:,:,xsx-plume.cs:xsx+plume.cs,:],2)
    #for 'smoke' create cross-section totals
    qvapor_mtot = np.nansum(qvapor[:,:,xsx-plume.cs:xsx+plume.cs,:],2)
