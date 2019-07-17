#inpue module for the original dry run for RxCADRE (initial submission to Atmosphere July 2019)
import numpy as np

#paths
wrfpath = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs_spinup'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile


#time setup of the domain
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)


#spatial setup of the domain
sfc_hgt = 62 							#surface height MSL (m)
ll_utm = np.array([517000,3377000])     #lower left corner of domain


#observational data
bounds_shape = '/Users/nmoisseeva/data/qgis/LG2012_WGS'

#ignition setup
fire_dict_utm = {'fireline1':{'start':np.array([525828,3379011]), 'end':np.array([524551,3378179])},\
				'fireline2':{'start':np.array([525729,3379075]), 'end':np.array([524487,3378275])},\
				'fireline3':{'start':np.array([525612,3379181]), 'end':np.array([524409,3378388])},\
				'fireline4':{'start':np.array([525549,3379284]), 'end':np.array([524331,3378480])} }
fuel_cat = 1
