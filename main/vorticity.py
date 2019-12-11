#December 2019
#nmoisseeva@eoas.ubc.ca



import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import sys
import imp
import pickle
from scipy.io import netcdf


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

#=================end of input===============


for nCase,Case in enumerate(plume.fireline_runs):

    interppath = plume.wrfdir + 'interp/wrfinterp_' + Case + '.npy'
    print('Opening data: %s' %interppath)
    # interpdict = np.load(interppath, allow_pickle=True).item()
    interpfile = open(interppath,'rb')
    interpdict = pickle.load(interpfile)   # load here the above pickle

    dimT, dimZ, dimY, dimX = np.shape(interpdict['U'])

    #dump netcdf file for ParaView
    netcdfout = netcdf.netcdf_file(plume.wrfdir + 'velField_%s.nc' %Case, 'w')
    netcdfout.history = 'Destaggered interpolated velocity field'

    dimTime = np.arange(0,dimT,4)

    netcdfout.createDimension('Time', len(dimTime))
    Time = netcdfout.createVariable('Time', 'i', ('Time',))
    Time[:] = np.arange(len(dimTime))
    Time.units = 'minutes'

    netcdfout.createDimension('dimz', dimZ)
    dimz = netcdfout.createVariable('dimz', 'i', ('dimz',))
    dimz[:] = plume.lvl
    dimz.units = 'meters'

    netcdfout.createDimension('dimy', dimY)
    dimy = netcdfout.createVariable('dimy', 'i', ('dimy',))
    dimy[:] = np.arange(0, dimY*plume.dy, plume.dy)
    dimy.units = 'meters'

    netcdfout.createDimension('dimx', int(dimX/3))
    dimx = netcdfout.createVariable('dimx', 'i', ('dimx',))
    dimx[:] = np.arange(0, int(dimX/3)*plume.dx, plume.dx)
    dimx.units = 'meters'

    u = netcdfout.createVariable('U', 'f', ('Time','dimz','dimy','dimx',))
    u.units = 'meters/sec'
    v = netcdfout.createVariable('V', 'f', ('Time','dimz','dimy','dimx',))
    v.units = 'meters/sec'
    w = netcdfout.createVariable('W', 'f', ('Time','dimz','dimy','dimx',))
    w.units = 'meters/sec'
    grnhfx = netcdfout.createVariable('grnhfx', 'f', ('Time','dimy','dimx',))
    grnhfx.units = 'Watt/meters^2'

    for nTime, time in enumerate(dimTime):
        u[nTime,:,:,:] = interpdict['U'][time,:,:,:int(dimX/3)]
        v[nTime,:,:,:]  = interpdict['V'][time,:,:,:int(dimX/3)]
        w[nTime,:,:,:] = interpdict['W'][time,:,:,:int(dimX/3)]
        grnhfx[nTime,:,:] = interpdict['GRNHFX'][time,:,:int(dimX/3)]
    netcdfout.close()
