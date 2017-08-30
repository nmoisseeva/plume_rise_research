#script to convert micromet data to sensible heat flux
#nmoisseeva@eoas.ubc.ca
#Aug 2017

import numpy as np

datapath = '/Users/nmoisseeva/Desktop/Htest.csv'
data = np.genfromtxt(datapath,usecols=(1,2,3,4),skip_header=4,delimiter=',',dtype=float)

fs = 0.170
fm = 0.02218

Ts_bar = np.mean(data[:,1])
Uz_bar = np.mean(data[:,0])
# Qt_bar = np.mean(data[:,2])/fs
# Qr_bar = np.mean(data[:,3])/fm
# Hs = Qt_bar - Qr_bar
# print Hs/(1005*1.2)


Ts_prime = data[:,1] - Ts_bar
Uz_prime = data[:,0] - Uz_bar

Cov = Ts_prime*Uz_prime

aveCov = np.mean(Cov)
H = aveCov*1005*1.2
print H

