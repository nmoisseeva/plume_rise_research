# nmoisseeva@eoas.ubc.ca
# June 2018


import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import sys
import imp
from matplotlib import animation


#====================INPUT===================

#all common variables are stored separately
import plume_q as plume
imp.reload(plume) 	#force load each time
#=================end of input===============


for nCase,Case in enumerate(plume.tag):

    cspath = plume.wrfdir + 'interp/wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    dimT, dimZ, dimX = np.shape(csdict['temp'])

    #save initial temperature prfile
    profpath = plume.wrfdir + 'interp/profT0' + Case + '.npy'
    profileT = np.mean(csdict['temp'][0,:,:],1)
    np.save(profpath,profileT)


    fig = plt.figure(figsize=(12,6))
    for nTime,Time in enumerate(np.arange(0,dimT,15)):
        plt.plot(csdict['w'][Time,1,:], label='%d min' %(Time*15/60.))
        plt.legend()
        plt.title('NEAR-SURFACE W: %s' %Case)
        ax = plt.gca()
        ax.set_xlim([0,dimX])
        ax.set_xticks(np.arange(0,dimX,int(dimX/10)))
        ax.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        plt.xlabel('horizontal distance [km]')
        plt.ylabel('vertical velocity w [m/s]')
    # plt.show()
    plt.savefig(plume.figdir + 'cs_analysis/surfW_%s.pdf' %Case)
    print('.....surface W plots saved in: %s' %(plume.figdir + 'cs_analysis/surfW_%s.pdf' %Case))
    plt.close()


    #create an animation of horizontal velocity---------------------------------------------
    print('.....creating vertical crossection of U + H2O animation')
    fig = plt.figure(figsize=(8,6))
    ax = plt.gca()
    ax.set_title('HORIZONTAL VELOCITY CONVERGENCE: %s' %Case)
    maxU = int(np.max(csdict['u']))
    maxQ = int(np.max(csdict['qvapor']))
    uLevels = np.arange(-maxU,maxU+.1,maxU/20.)
    qLevels = np.arange(maxQ*0.05,maxQ,maxQ/40.)
    # create initial frame
    # ---u contours and colorbar
    cntrf = ax.contourf(csdict['u'][0,:,:], cmap=plt.cm.Spectral_r,levels=uLevels,extend='both')
    cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
    cbarf.set_label('horizontal velocity $[m s^{-1}]$')
    ax.set_xlabel('horizontal distance [km]')
    ax.set_ylabel('height AGL [m]')
    ax.set_xlim([0,dimX])
    # ---non-filled vapor contours and colorbar
    cntr = ax.contour(csdict['qvapor'][0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
    #
    # ---heat flux
    axh = ax.twinx()
    axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh.set_ylim([0,140])
    axh.set_xlim([0,dimX])
    axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
    axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
    axh.tick_params(axis='y', colors='red')
    ln = axh.plot(csdict['ghfx'][0,:], 'r-')
    # fig.tight_layout()


    def update_plot(n,csdict,cntrf,cntr):
        ax.clear()
        ax.set_title('HORIZONTAL VELOCITY CONVERGENCE: %s' %Case)
        cntrf = ax.contourf(csdict['u'][n,:,:],cmap=plt.cm.Spectral_r, levels=uLevels,extend='both')
        cntr = ax.contour(csdict['qvapor'][n,:,:], cmap=plt.cm.Greys,levels=qLevels,linewidths=0.6)
        ax.set_xlabel('horizontal distance [km]')
        ax.set_ylabel('height AGL [m]')
        ax.set_yticks(np.arange(0,len(plume.lvl),10))
        ax.set_yticklabels(plume.lvl[::10])
        axh.clear()
        axh.set_ylim([0,140])
        axh.set_xlim([0,dimX])
        axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        ln = axh.plot(csdict['ghfx'][n,:], 'r-')
        axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        # fig   .tight_layout()

        return cntrf, ln, cntr,

    #plot all frames
    ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,cntrf,cntr), interval=3)
    # plt.show()
    ani.save(plume.figdir + 'anim/u/u%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
    plt.close()
    print('.....saved in: %s' %(plume.figdir + 'anim/u/u%s.mp4' %Case))

    #create an animation of vertical velocity---------------------------------------------
    print('.....creating vertical crossection of W + H2O animation')
    fig = plt.figure(figsize=(8,6))
    ax = plt.gca()
    ax.set_title('VERTICAL VELOCITY : %s' %Case)
    maxW = int(np.max(csdict['w']))
    maxQ = int(np.max(csdict['qvapor']))
    wLevels = np.arange(-maxW,maxW+.1,maxW/20)
    qLevels = np.arange(maxQ*0.05,maxQ,maxQ/40)
    # create initial frame
    # ---u contours and colorbar
    cntrf = ax.contourf(csdict['w'][0,:,:], cmap=plt.cm.PRGn_r,levels=wLevels,extend='both')
    cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
    cbarf.set_label('vertical velocity $[m s^{-1}]$')
    ax.set_xlabel('horizontal distance [km]')
    ax.set_ylabel('height AGL [m]')
    ax.set_xlim([0,dimX])
    # ---non-filled vapor contours and colorbar
    cntr = ax.contour(csdict['qvapor'][0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
    #
    # ---heat flux
    axh = ax.twinx()
    axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    axh.set_ylim([0,140])
    axh.set_xlim([0,dimX])
    axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
    axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
    axh.tick_params(axis='y', colors='red')
    ln = axh.plot(csdict['ghfx'][0,:], 'r-')
    # fig.tight_layout()


    def update_plot(n,csdict,cntrf,cntr):
        ax.clear()
        ax.set_title('VERTICAL VELOCITY: %s' %Case)
        cntrf = ax.contourf(csdict['w'][n,:,:],cmap=plt.cm.PRGn_r, levels=wLevels,extend='both')
        cntr = ax.contour(csdict['qvapor'][n,:,:], cmap=plt.cm.Greys,levels=qLevels,linewidths=0.6)
        ax.set_xlabel('horizontal distance [km]')
        ax.set_ylabel('height AGL [m]')
        ax.set_yticks(np.arange(0,len(plume.lvl),10))
        ax.set_yticklabels(plume.lvl[::10])
        axh.clear()
        axh.set_ylim([0,140])
        axh.set_xlim([0,dimX])
        axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        ln = axh.plot(csdict['ghfx'][n,:], 'r-')
        axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        # fig   .tight_layout()
        return cntrf, ln, cntr,

    #plot all frames
    ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,cntrf,cntr), interval=3)
    # plt.show()
    ani.save(plume.figdir + 'anim/w/w%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
    plt.close()
    print('.....saved in: %s' %(plume.figdir + 'anim/w/w%s.mp4' %Case))
