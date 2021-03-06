#inpue module for the original dry run for RxCADRE (initial submission to Atmosphere July 2019)
import numpy as np
#paths
wrfdata = '/Users/nmoisseeva/data/plume/rxcadre/wrfout_L2G_cat1_obs'
fig_dir = '/Users/nmoisseeva/code/plume/figs/RxCADRE/'
interp_path = '/Users/nmoisseeva/code/plume/RxCADRE/npy/qv_L2G_cat1obs_interp.npy'
spinup_path = '/Users/nmoisseeva/data/plume/RxCADRE/Feb2019/wrfout_L2G_cat1obs_spinup'
# spinup_path = '/Users/nmoisseeva/data/plume/RxCADRE/wrfout_L2G_cat1obs_spinup'
wrfinput = '/Users/nmoisseeva/data/plume/rxcadre/wrfinput_d01'
input_fc = '/Users/nmoisseeva/sfire/wrf-fire/WRFV3/test/em_fire/rxcadre_moist/input_fc'


#time setup of the domain
runstart = '10:00:00' 					#start time (if restart run time of inital simulation)
run_min = 45 							#run length in minutes
hist_int = 10 							#history interval in seconds (main run)
spinup_hist_int = 60                   #spinup history interval (seconds)
moist_run = 0 							#no moisture included (old simulations)

#spatial setup of the domain
sfc_hgt = 62 							#surface height MSL (m)
ll_utm = np.array([517000,3377000])     #lower left corner of domain
basemap_path = '/Users/nmoisseeva/code/plume/rxcadre/npy/%s_%s_bm_fire.npy' %(ll_utm[0],ll_utm[1])
lvl = np.arange(0,1700,20)

#observational data
bounds_shape = '/Users/nmoisseeva/data/rxcadre/qgis/BurnPerimeterWGS'
disp_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/SmokeDispersion_L2G_20121110.csv'
emis_data = '/Users/nmoisseeva/data/RxCADRE/dispersion/Data/Emissions_L2G_20121110.csv'
pre_moisture = '/Users/nmoisseeva/data/RxCADRE/meteorology/soundings/MoistureProfile_NM.csv'            #pre-burn moisture profile



#ignition setup (LWIR estimate)
fire_dict_utm = {'fireline1':{'start':np.array([525828,3379011]), 'end':np.array([524551,3378179])},\
				'fireline2':{'start':np.array([525729,3379075]), 'end':np.array([524487,3378275])},\
				'fireline3':{'start':np.array([525612,3379181]), 'end':np.array([524409,3378388])},\
				'fireline4':{'start':np.array([525549,3379284]), 'end':np.array([524331,3378480])} }
fuel_cat = 1


#instrument locations (UTM)
csu_lcn = [525803.12, 3378544.15]
met_lcn = [525128.85, 3378208.91]
rad_lcn = [526090,3378766]              #MUST DOUBLE-CHECK CONVERSION TO UTM


#flight profiles
profiles_start = ['12:02:30','12:13:00','12:27:00']
profiles_end = ['12:05:30','12:19:00','12:29:30']
garage = ['12:36:00','13:04:00']
corkscrew = ['13:06:30','13:12:00']
