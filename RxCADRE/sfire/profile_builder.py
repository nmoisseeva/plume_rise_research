#profile builder - interpolates and averages two soundings

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


profile_data = '/Users/nadya2/code/plume/RxCADRE/sfire/two_profiles_converted.csv'
save_data = '/Users/nadya2/code/plume/RxCADRE/sfire/input_sounding'
lvl = np.arange(0,4000,10)

#import profiles
print('Importing dispersion data from %s' %profile_data)
disp_dict = {}
am = np.genfromtxt(profile_data, skip_header=4, usecols = [0,1,6,7], delimiter=',')
pm = np.genfromtxt(profile_data, skip_header=4, usecols = [10,11,16,17], delimiter=',')
# z10,T10,u10,v10 = np.genfromtxt(profile_data, skip_header=4, usecols = [0,1,6,7], delimiter=',')
# z2,T2,u2,v2 = np.genfromtxt(profile_data, skip_header=4, usecols = [10,11,16,17], delimiter=',')


fT10 = interpolate.interp1d(am[:,0],am[:,1],fill_value="extrapolate")
fu10 = interpolate.interp1d(am[:,0],am[:,2],fill_value="extrapolate")
fv10 = interpolate.interp1d(am[:,0],am[:,3],fill_value="extrapolate")
fT2 = interpolate.interp1d(pm[:,0],pm[:,1],fill_value="extrapolate")
fu2 = interpolate.interp1d(pm[:,0],pm[:,2],fill_value="extrapolate")
fv2 = interpolate.interp1d(pm[:,0],pm[:,3],fill_value="extrapolate")

iT10,iu10,iv10 = fT10(lvl),fu10(lvl),fv10(lvl)
iT2,iu2,iv2 = fT2(lvl),fu2(lvl),fv2(lvl)

ave_sounding = np.empty((len(lvl),5))
ave_sounding[:,0] = lvl
ave_sounding[:,1] = (iT10+iT2)/2.
ave_sounding[:,2] = 0.
ave_sounding[:,3] = (iu10+iu2)/2.
ave_sounding[:,4] = (iv10+iv2)/2.

np.savetxt(save_data, ave_sounding, fmt='%.2f', delimiter=' ', newline='\n')