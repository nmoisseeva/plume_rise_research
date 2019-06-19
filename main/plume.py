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

cs = 20                         #+/- grids for cross-section
wi, wf = 10, 65


# dirlist = os.listdir(wrfdir+'interp/') 	#get all files in directory
dirpath = wrfdir+'interp/wrfave_*'
dirlist = glob.glob(dirpath) #get all  interp files in directory
tag = [i[len(dirpath)-1:-4] for i in dirlist]    #W*S*F*R0
# tag = ['W8S400F7R0']

#exclude list (low wind speeds, F1, F4 (too hot), too high, randomly low, not sure)
exclude_runs = ['W1S400F7R0','W2S400F3R0','W2S400F1R0',\
                'W4S400F1R0','W12S400F1R0','W8S400F1R0',\
                'W7S0F4R0','W3S0F4R0','W5S0F4R0','W9S0F4R0',\
                'W4S400F13R0',\
                'W4S400F9R0','W4S400F8R0','W4Sn200F8R0',\
                'W11S400F10R0']


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
        elif str_tag=='S':
            if 'Sn' in letters:
                out_array.append(int(numbers[1])*-1)
            else:
                out_array.append(int(numbers[1]))
        elif str_tag=='F':
            out_array.append(int(numbers[2]))
        elif str_tag=='R':
            out_array.append(int(numbers[3]))
    out_array = np.array(out_array)
    return out_array
