# January 2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time


#=================end of input===============

RunList = [i for i in plume.tag if i not in plume.exclude_runs]
runCnt = len(RunList)


#take a vertical cross-section (ensure it's stable)
#calculate CDF of cumT (delT*dz)
#calculate temperature gradient -
#map temperature gradient to CDF
