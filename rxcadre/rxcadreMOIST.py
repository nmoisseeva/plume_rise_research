#inpue module for the omoist RxCADRE run
import numpy as np

#paths
wrfdata = '/Users/nmoisseeva/data/plume/rxcadre/wrfout_test_main'
fig_dir = '/Users/nmoisseeva/code/plume/rxcadre/figs/'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'
geo_path = '/Users/nmoisseeva/code/plume/rxcadre/npy/wrfout_geo.npy'
spinup_path = '/Users/nmoisseeva/data/plume/rxcadre/wrfout_sfmo50_spinup'
wrfinput = '/Users/nmoisseeva/data/plume/rxcadre/wrfinput_d01'
input_fc = '/Users/nmoisseeva/sfire/wrf-fire/WRFV3/test/em_fire/rxcadre_moist/input_fc'


#time setup of the domain
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)
run_min = 45 							#run length in minutes
hist_int = 10 							#history interval in seconds (main run)
spinup_hist_int = 60                   #spinup history interval (seconds)
moist_run = 1 							#moisture included (new simulations)


#spatial setup of the domain
sfc_hgt = 62 							#surface height MSL (m)
ll_utm = np.array([517000,3377000])     #lower left corner of domain
basemap_path = '/Users/nmoisseeva/code/plume/rxcadre/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])

#observational data
bounds_shape = '/Users/nmoisseeva/code/plume/rxcadre/gis/BurnPerimeterWGS'
hip1_shape = '/Users/nmoisseeva/code/plume/rxcadre/gis/HIP1SensorLocationsWGS'
# hip1_shape = '/Users/nmoisseeva/data/rxcadre/instruments/RDS-2016-0014/Data/SurveyPointInformation'
disp_data = '/Users/nmoisseeva/data/rxcadre/dispersion/RDS-2014-0015/Data/SmokeDispersion_L2G_20121110.csv'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv' #pre-burn moisture profile


#ignition setup - GPS based
fire_dict_utm = {'fireline1':{'start':np.array([525507,3379248]), 'end':np.array([524305,3378463])},\
				'fireline2':{'start':np.array([525630,3379179]), 'end':np.array([524384,3378366])},\
				'fireline3':{'start':np.array([525732,3379077]), 'end':np.array([524462,3378248])},\
				'fireline4':{'start':np.array([525820,3379004]), 'end':np.array([524548,3378172])} }
fuel_cat = 1


#instrument locations (UTM)
# csu_lcn = [525803.12, 3378544.15]         #don't think data is from here, but from the met tower
met_lcn = [525128.85, 3378208.91]           #20-ft micromet tower ID7001, FLag 1502 (within lot)
# rad_lcn = [526090,3378766]              #MUST DOUBLE-CHECK CONVERSION TO UTM

#flight profiles
profiles_start = ['12:02:30','12:13:00','12:27:00']
profiles_end = ['12:05:30','12:19:00','12:29:30']
garage = ['12:36:00','13:04:00']
corkscrew = ['13:06:30','13:12:00']
