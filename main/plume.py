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
dy = 40.

cs = 10                         #+/- grids for cross-section
wi, wf = 25, 375
fireline = 4000.,6000.          #fireline start and end in meters
ign_over = 20                   #number of history intervals exluded from start

# dirlist = os.listdir(wrfdir+'interp/') 	#get all files in directory
dirpath = wrfdir+'interp/wrfave_*'
dirlist = glob.glob(dirpath) #get all  interp files in directory
# tag = [i[len(dirpath)-1:-4] for i in dirlist]    #W*S*F*R0
tag = ['W4F7R4']

#exclude list (F1, F8?,F9?)
exclude_runs = ['W5F4R0','W5F4R1','W5F4R2','W5F4R3','W5F4R4' ]
# fireline_runs = ['W4F7R4L1']
fireline_runs = ['W4F7R4L1','W4F7R4','W4F7R4L4']


#common functions
#----------------------------------------------------------------------------
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
            if letters[-2]=='L':        #-2 accounts for empty string at the end of the file
                out_array.append(int(numbers[3]))
            else:
                out_array.append(2)
    out_array = np.array(out_array)
    return out_array

#load crosssection dictionary, extract profiles them and load them in to variables
def load_CS_prep_Profiles(Case):
    import numpy as np
    #load cross section
    cspath = wrfdir + 'interp/wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    #save initial profiles
    profpathT = wrfdir + 'interp/profT0' + Case + '.npy'
    profpathU = wrfdir + 'interp/profU0' + Case + '.npy'

    if not os.path.isfile(profpathT) or not os.path.isfile(profpathU):
        profileT = np.mean(csdict['temp'][0,:,:],1)
        np.save(profpathT,profileT)

        profileU = np.mean(csdict['u'][0,:,:],1)
        np.save(profpathU,profileU)
        print('...Generated new profiles' )
    return csdict

#get BL top
def get_zi(T0):
    dT = T0[1:]-T0[0:-1]
    gradT = dT[1:] - dT[0:-1]
    si = 3
    zi_idx = np.argmax(gradT[si:]) + si                 #vertical level index of BL top
    zi = dz * zi_idx
    return zi
