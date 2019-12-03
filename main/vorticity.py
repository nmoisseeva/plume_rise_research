#December 2019
#nmoisseeva@eoas.ubc.ca



import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import sys
import imp
from matplotlib import animation
from numpy import ma
from matplotlib import ticker
import metpy


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

doAnimations = 0    #flag to do animations:
#=================end of input===============


for nCase,Case in enumerate(plume.fireline_runs):

    interppath = plume.wrfdir + 'interp/wrfinterp_' + Case + '.npy'
    print('Opening data: %s' %interppath)
    interpdict = np.load(interppath, allow_pickle=True).item()

    dimT, dimZ, dimY, dimX = np.shape(interpdict['u'])

    vorticity = mepty.calc.vorticity(interpdict['u'][70, 1,:,:], interpdict['v'][70, 1,:,:], plume.dx, plume.dy)
    plt.contourf(vorticity)
    plt.colorbar()
    plt.show()


    if doAnimations:

        vorticityArray = np.empty((dimT,dimY,dimX))
        for nT in dimT:
            vorticityArray[nT,:,:] = mepty.calc.vorticity(interpdict['u'][nT,1,:,:], interpdict['v'][nT, 1,:,:], plume.dx, plume.dy)


        #create an animation of vorticity---------------------------------------------
        fig = plt.figure(figsize=(10,5))
        gs = fig.add_gridspec(ncols=2, nrows=1,width_ratios=[3,1])
        print('.....creating animation of surface vertical vorticity')
        plt.suptitle('VERTICAL VORTICITY: %s' %Case)

        ax1=fig.add_subplot(gs[0])
        # create initial frame
        # ---u contours and colorbar
        im = ax1.imshow(vorticityArray[0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.dy],cmap=plt.cm.RdBu_r)
        ax1.set_aspect('auto')
        cbari = fig.colorbar(im, orientation='horizontal')
        cbari.set_label('vertical vorticity $[s^{-1}]$')

        ax2=fig.add_subplot(gs[1])
        fim = ax2.imshow(interpdict['ghfx'][0,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
        ax2.set_aspect('auto')
        cbarif = fig.colorbar(fim, orientation='horizontal')
        cbarif.set_label('heat flux [$kW / m^2$]')
        ax2.set_xlabel('x distance [m]')
        ax2.set_ylabel('y distance [m]')
        # plt.subplots_adjust(top=0.85)
        plt.tight_layout(rect=[0, 0, 1, 0.92])

        def update_plot(n,csdict,cntrf,cntr):
            ax1.clear()
            im = ax1.imshow(vorticityArray[n,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.dy],cmap=plt.cm.RdBu_r)
            ax1.set_aspect('auto')
            ax1.set_xlabel('x [m]')
            ax1.set_ylabel('y[m]')
            ax1.set_ylim([0,plume.dimY])
            axh1.clear()
            axh1.set_ylim([0,dimY*plume.dy])
            axh1.set_xlim([0,dimX*plume.dx])

            ax2.clear()
            fim = ax2.imshow(interpdict['ghfx'][n,75:175,0:75],cmap=plt.cm.YlOrRd,extent=[0,3000,3000,7000],vmin=0, vmax = 150)
            ax2.set_aspect('auto')
            ax2.set_xlabel('x distance [m]')
            ax2.set_ylabel('y distance [m]')

            return cntrf, ln, cntr,

        #plot all frames
        ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,im,cntr), interval=3)
        # plt.show()
        ani.save(plume.figdir + 'anim/vorticity/curl%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        plt.close()
        print('.....saved in: %s' %(plume.figdir + 'anim/vorticity/curl%s.mp4' %Case))
