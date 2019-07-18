#inpue module for the omoist RxCADRE run
import numpy as np

#paths
wrfdata = '/Users/nmoisseeva/data/plume/rxcadre/wrfout_L2G_cat1obs'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'


#time setup of the domain
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)
run_min = 45 							#run length in minutes
hist_int = 10 							#history interval in seconds


#spatial setup of the domain
sfc_hgt = 62 							#surface height MSL (m)
ll_utm = np.array([517000,3377000])     #lower left corner of domain

#observational data
bounds_shape = '/Users/nmoisseeva/data/qgis/LG2012_WGS'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile



#ignition setup - GPS based
fire_dict_utm = {'fireline1':{'start':np.array([525507,3379248]), 'end':np.array([524305,3378463])},\
				'fireline2':{'start':np.array([525630,3379179]), 'end':np.array([524384,3378366])},\
				'fireline3':{'start':np.array([525732,3379077]), 'end':np.array([524462,3378248])},\
				'fireline4':{'start':np.array([525820,3379004]), 'end':np.array([524548,3378172])} }
fuel_cat = 1


#instrument locations (UTM)
csu_lcn = [525803.12, 3378544.15]
met_lcn = [525128.85, 3378208.91]
