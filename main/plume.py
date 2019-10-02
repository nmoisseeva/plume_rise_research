#nmoisseeva@eoas.ubc.ca
#March 2019

import numpy as np
import os.path
import glob

#input values for plume analysis
#--------------------------------
wrfdir = '/Users/nmoisseeva/data/plume/main/'
figdir = '/Users/nmoisseeva/code/plume/main/figs/'
dz = 40
lvl = np.arange(0,2800,dz)	 	#vertical levels in m
dx = 40.                        #horizontal grid spacing

cs = 10                         #+/- grids for cross-section
wi, wf = 10, 65


# dirlist = os.listdir(wrfdir+'interp/') 	#get all files in directory
dirpath = wrfdir+'interp/wrfave_*'
dirlist = glob.glob(dirpath) #get all  interp files in directory
tag = [i[len(dirpath)-1:-4] for i in dirlist]    #W*S*F*R0
# tag = ['W8S400F7R0']

#exclude list (low wind speeds, F1, F4 (too hot), too high, randomly low, not sure)
exclude_runs = []
fireline_runs = ['W5F7R2','W5F7R2L2','W5F7R2L3','W5F7R2L410km']

#common functions
#--------------------------------
#tag reading function read_tag(variable type, string array)
def read_tag(str_tag, str_array):
    import re
    out_array = []
    for nTag,tag in enumerate(str_array):
        letters = re.split('\d+', tag)
        numbers = re.findall('\d+', tag)
        if str_tag=='W':
            out_array.append(int(numbers[0]))
        elif str_tag=='F':
            out_array.append(int(numbers[1]))
        elif str_tag=='R':
            out_array.append(int(numbers[2]))
        elif str_tag=='L':
            out_array.append(int(numbers[3]))
    out_array = np.array(out_array)
    return out_array
