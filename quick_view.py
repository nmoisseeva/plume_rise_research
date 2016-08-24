from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib as mpl 
import warnings
warnings.filterwarnings("ignore")

# wrfdata='/Users/nadya2/data/plume/comps/D/wrfout_b'
# wrfdata='/Users/nadya2/data/plume/sensitivity/wind/wrfout_6v'


nc_data = netcdf.netcdf_file(wrfdata, mode ='r')   
