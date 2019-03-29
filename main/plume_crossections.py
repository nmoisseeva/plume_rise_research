#nmoisseeva@eoas.ubc.ca
#March 2019

import numpy as np
from scipy.io import netcdf
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
    wrfpath = plume.wrfdir + 'wrfout_'+ Case
    wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')
    ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX'))
    wrfdata.close()

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
    u = interpdict['U']

    #get dimensions
    dimt, dimy, dimx = np.shape(ghfx)
    xsx = int(round(dimy/2.))

    var_list = ['ghfx','qvapor','temp','w','u']
    csdict = {}

    for variable in var_list:
        print(variable)
        #create fire cross-section
        if variable == 'ghfx':
            fire_ave = np.mean(vars()[variable][:,:,xsx-plume.cs:xsx+plume.cs,:],1)
            csdict[variable] = fire_ave
            xmax = np.argmax(fire_ave,axis=1)
        elif variable == 'qvapor':
            fire_tot = np.nansum(vars()[variable][:,:,xsx-plume.cs:xsx+plume.cs,:],2)
            csdict[variable] = fire_tot
        else:
            fire_ave = np.mean(vars()[variable][:,:,xsx-plume.cs:xsx+plume.cs,:],2)
            csdict[variable] = fire_ave

    #create time-average around peak flux--------------------------
    tdict = {}

    for variable in var_list:
        tvar = []
        for nP, pt in enumerate(xmax[1:]):
            tvar.append(csdict[variable][nP+1,pt-plume.wi:pt+plume.wf])
        tdict[variable] = tvar
