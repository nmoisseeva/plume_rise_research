#nmoisseeva@eoas.ubc.ca
#March 2019

import numpy as np
import os.path

#input values for plume analysis
#--------------------------------
wrfdir = '/Users/nmoisseeva/data/plume/main/'
figdir = '/Users/nmoisseeva/code/plume/main/figs/'
lvl = np.arange(0,2500,40)	 	#vertical levels in m
dx = 40.                        #horizontal grid spacing

cs = 20                         #+/- grids for cross-section                               

dirlist = os.listdir(wrfdir+'interp/') 	#get all files in directory
tag = [i[7:-4] for i in dirlist[1:]]    #W*S*F*R0

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
    out_array = np.array(out_array)
    return out_array
