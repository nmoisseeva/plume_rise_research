from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib as mpl 
import warnings
warnings.filterwarnings("ignore")

# wrfdata='/Users/nadya2/data/plume/comps/D/wrfout_b'
wrfdata='/Users/nadya2/data/plume/sensitivity/wind/wrfout_6v'
sounding = '/Users/nadya2/code/plume/comps/B/input_sounding_a'

nc_data = netcdf.netcdf_file(wrfdata, mode ='r')   


x = np.genfromtxt(sounding, dtype=float, skip_header=1,usecols=(0,1))
plt.figure(figsize=(6,6))
plt.plot(x[:,1], x[:,0])
plt.xlim([295,315])
plt.ylim([0,2000])
plt.ylabel('height [m]')
plt.xlabel('potential temperature [K]')
plt.tight_layout()
plt.savefig('/Users/nadya2/GoogleDrive/PhD/comps_written/roland/figs/sounding.pdf')