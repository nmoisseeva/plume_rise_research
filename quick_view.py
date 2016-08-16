from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib as mpl 
import warnings
warnings.filterwarnings("ignore")

wrfdata='/Users/nadya2/data/plume/comps/B/wrfout_a'

nc_data = netcdf.netcdf_file(wrfdata, mode ='r')   
