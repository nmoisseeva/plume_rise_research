# nmoisseeva@eoas.ubc.ca
# June 2018


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import sys
import imp
from matplotlib import animation
from numpy import ma
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#====================INPUT===================

#all common variables are stored separately
import Plume
imp.reload(Plume) 	#force load each time
import utils
import config

#=================end of input===============
Case = 'W4F7R4L4'


csdict = utils.prepCS(Case)
T0 = np.load(config.wrfdir + 'profiles/profT0' + Case + '.npy')
U0 = np.load(config.wrfdir + 'profiles/profU0' + Case + '.npy')



dimT, dimZ, dimX = np.shape(csdict['pm25'])
dimY = np.shape(csdict['ghfx2D'])[1]


# maxPM = int(np.max(csdict['pm25']))
# PMcontours = ma.masked_where(csdict['pm25'] <= 30,csdict['pm25'])


pmlvls = np.arange(0,dimZ*config.dz,config.dz)
cropX = int(dimX*0.7)
axMax = cropX * config.dx
haxis = np.arange(cropX)*config.dx
PMmg = csdict['pm25']/1000.                                        #smoke concentration in ppm
maxPM = int(np.max(PMmg))

PMcontours = ma.masked_where(PMmg <= 0.03,PMmg)
pmLevels = np.geomspace(0.03,maxPM/(10),num=10)


fig = plt.figure(figsize=(13,8))
gs = fig.add_gridspec(ncols=2, nrows=3,width_ratios=[7,1])

# create initial frame
ax1=fig.add_subplot(gs[0,:])
# ---cwi smoke  and colorbar
im = ax1.imshow(PMmg[0,:,:cropX], origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.cubehelix_r,vmin=0, vmax=maxPM/10)
cbaxes = inset_axes(ax1, width="30%", height="5%", loc=1)
cbari = fig.colorbar(im, cax=cbaxes, orientation='horizontal',label='CWI smoke [mg/kg]')
ax1.set(ylabel='height [m]', aspect='equal')

# ax2=fig.add_subplot(gs[0,-1])
# fim = ax2.imshow(csdict['ghfx2D'][0,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
# cbaxesf = inset_axes(ax2, width="40%", height="5%", loc=1)
# cbarif = fig.colorbar(fim, cax=cbaxesf, orientation='horizontal',label='$kW / m^2$')
# ax2.set(xlabel = 'x distance [m]', ylabel='y distance [m]',aspect='equal')

ax3=fig.add_subplot(gs[1,:])
uim = ax3.imshow((csdict['u'][0,:,:cropX].T - U0).T, origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.RdBu_r,vmin=-6, vmax=6)
ucntr = ax3.contour(PMcontours[0,:,:cropX],extent=[0,axMax,0,pmlvls[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
cbaxesu = inset_axes(ax3, width="30%", height="5%", loc=1)
cbariu = fig.colorbar(uim, cax=cbaxesu, orientation='horizontal',label='ralative horizontal velocity [m/s]')
ax3.set( ylabel='height [m]',aspect='equal')


ax4=fig.add_subplot(gs[2,:])
wim = ax4.imshow(csdict['w'][0,:,:cropX], origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.PRGn_r, vmin=-4, vmax=4)
wcntr = ax4.contour(PMcontours[0,:,:cropX],extent=[0,axMax,0,pmlvls[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
cbaxesu = inset_axes(ax4, width="30%", height="5%", loc=1)
cbariu = fig.colorbar(wim, cax=cbaxesu, orientation='horizontal',label='relative vertical velocity [m/s]')
ax4.set(xlabel = 'x distance [m]', ylabel='height [m]',aspect='equal')

def update_plot(n,PMmg,im,uim,ucntr,wim,wcntr):
    ax1.clear()
    im = ax1.imshow(PMmg[n,:,:cropX], origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.cubehelix_r,vmin=0, vmax=maxPM/10)
    ax1.set(ylabel='height [m]',aspect='equal')

    # ax2.clear()
    # fim = ax2.imshow(csdict['ghfx2D'][n,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,2500,7500],vmin=0, vmax = 150)
    # ax2.set(xlabel = 'x distance [m]', ylabel='y distance [m]',aspect='equal')

    ax3.clear()
    uim = ax3.imshow((csdict['u'][n,:,:cropX].T - U0).T, origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.RdBu_r,vmin=-6, vmax=6)
    ucntr = ax3.contour(PMcontours[n,:,:cropX],extent=[0,axMax,0,pmlvls[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
    ax3.set(ylabel='height [m]',aspect='equal')

    ax4.clear()
    wim = ax4.imshow(csdict['w'][n,:,:cropX], origin='lower', extent=[0,axMax,0,pmlvls[-1]],cmap=plt.cm.PRGn_r, vmin=-4, vmax=4)
    wcntr = ax4.contour(PMcontours[n,:,:cropX],extent=[0,axMax,0,pmlvls[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
    ax4.set(xlabel = 'distance [m]',ylabel='height [m]',aspect='equal')

    return im, uim, ucntr, wim, wcntr,


#plot all frames
# ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,im,cntr), interval=3)
ani=animation.FuncAnimation(fig, update_plot, dimT,fargs=(PMmg,im,uim,ucntr,wim,wcntr), interval=3)


# ani.save(config.figdir + 'anim/w/smoke%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
ani.save(config.figdir + 'animation%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
plt.close()
