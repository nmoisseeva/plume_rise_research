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


#====================INPUT===================

#all common variables are stored separately
import plume
imp.reload(plume) 	#force load each time

doAnimations = 1    #flag to do animations:
#=================end of input===============


for nCase,Case in enumerate(plume.tag):

    cspath = plume.wrfdir + 'interp/wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    dimT, dimZ, dimX = np.shape(csdict['temp'])

    #save initial temperature prfile
    profpathT = plume.wrfdir + 'interp/profT0' + Case + '.npy'
    profileT = np.mean(csdict['temp'][0,:,:],1)
    np.save(profpathT,profileT)

    #save initial temperature prfile
    profpathU = plume.wrfdir + 'interp/profU0' + Case + '.npy'
    profileU = np.mean(csdict['u'][0,:,:],1)
    np.save(profpathU,profileU)

    # fig = plt.figure(figsize=(12,6))
    # for nTime,Time in enumerate(np.arange(0,dimT,15)):
    #     plt.plot(csdict['w'][Time,1,:], label='%d min' %(Time*15/60.))
    #     plt.legend()
    #     plt.title('NEAR-SURFACE W: %s' %Case)
    #     ax = plt.gca()
    #     ax.set_xlim([0,dimX])
    #     ax.set_xticks(np.arange(0,dimX,int(dimX/10)))
    #     ax.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
    #     plt.xlabel('horizontal distance [km]')
    #     plt.ylabel('vertical velocity w [m/s]')
    # # plt.show()
    # plt.savefig(plume.figdir + 'cs_analysis/surfW_%s.pdf' %Case)
    # print('.....surface W plots saved in: %s' %(plume.figdir + 'cs_analysis/surfW_%s.pdf' %Case))
    # # plt.close()
    #
    #
    #
    # #plot contours
    PMcontours = ma.masked_where(csdict['pm25'] <= 30,csdict['pm25'] )
    # fig = plt.figure(figsize=(18,10))
    # plt.suptitle('%s' %Case)
    # gs = fig.add_gridspec(ncols=3, nrows=2,height_ratios=[3,1])
    #
    # ax1=fig.add_subplot(gs[0,0])
    # # plt.subplot(2,3,1,gs[0])
    # plt.title('Time-averaged W')
    # # ---w contours and colorbar
    # im = ax1.imshow(avedict['w'], origin='lower',extent=[0,haxis[-1],0,plume.lvl[-1]], cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
    # ax1.set_aspect('auto')
    # cbari =fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.1)
    # cbari.set_label('horizontal velocity w $[m s^{-2}]$')
    # ax1.set_xlabel('horizontal distance [m]')
    # ax1.set_ylabel('height AGL [m]')
    # ax1.axhline(y = si * plume.dz, ls=':', c='lightgrey', label='surface layer height at ignition')
    # ax1.axhline(y = BLdict['zi'][nCase], ls=':', c='darkgrey', label='BL height at ignition')
    # ax1.axhline(y=plume_inj[nCase],ls='--', c='black',label='derived plume top')
    # ax1.axvline(x = sliceX*plume.dx, ls=':',c='black',label='location of concentration profile')
    # ax1.axvline(x = xmax*plume.dx, ls=':',c='darkgrey',label='tilt defintion')
    # ax1.legend()
    # # ---non-filled pm contours and colorbar
    # cntr = ax1.contour(PMcontours, extent=[0,haxis[-1],0,plume.lvl[-1]],locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
    # ax1.plot(haxis,centerline,ls='--', c='darkgrey',label='plume centerline' )
    # # ---heat flux
    # axh1 = ax1.twinx()
    # axh1.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
    # axh1.set_ylim([0,150])
    # axh1.tick_params(axis='y', colors='red')
    # ln = axh1.plot(haxis, avedict['ghfx'], 'r-')
    # axh1.set_xlim([0,haxis[-1]])
    #
    #
    #
    # ax2=fig.add_subplot(gs[0,1])
    # plt.title('Time-averaged U')
    # im = ax2.imshow((avedict['u'].T-U0[nCase]).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]], cmap=plt.cm.RdBu_r, vmin=-4, vmax=4)
    # ax2.set_aspect('auto')
    # # ---non-filled vapor contours and colorbar
    # cntr = ax2.contour(PMcontours, extent=[0,dimX*plume.dx,0,plume.lvl[-1]], locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
    # ax2.plot(np.arange(dimX)*plume.dx,centerline,ls='--', c='darkgrey' )
    # # ---heat flux
    # ln = axh2.plot(np.arange(dimX)*plume.dx, avedict['ghfx'], 'r-')
    # axh2.set_xlim([0,dimX*plume.dx])
    #
    #



    if doAnimations:
        #create an animation of horizontal velocity---------------------------------------------
        print('.....creating vertical crossection of U + PM2.5 animation')
        fig = plt.figure(figsize=(8,6))
        ax = plt.gca()
        ax.set_title('RELATIVE HORIZONTAL VELOCITY CONVERGENCE: %s' %Case)
        # create initial frame
        # ---u contours and colorbar
        im = ax.imshow((csdict['u'][0,:,:].T - profileU).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)

        ax.set_aspect('auto')
        cbari = fig.colorbar(im, orientation='horizontal',fraction=0.046, pad=0.1)
        cbari.set_label('ralative horizontal velocity $[m s^{-1}]$')
        # ---non-filled vapor contours and colorbar
        cntr = ax.contour(PMcontours[0,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
        ax.set_xlim([0,dimX*plume.dx])
        ax.set_xlabel('horizontal distance [km]')
        ax.set_ylabel('height AGL [m]')
        # ---heat flux
        axh = ax.twinx()
        ln = axh.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][0,:], 'r-')
        axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        axh.set_ylim([0,150])
        axh.set_xlim([0,dimX*plume.dx])
        axh.tick_params(axis='y', colors='red')


        def update_plot(n,csdict,cntrf,cntr):
            ax.clear()
            ax.set_title('RELATIVE HORIZONTAL VELOCITY CONVERGENCE: %s' %Case)
            im = ax.imshow((csdict['u'][n,:,:].T - profileU).T,origin='lower',extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)
            cntr = ax.contour(PMcontours[n,:,:], extent=[0,dimX*plume.dx,0,plume.lvl[-1]],locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=1)
            ax.set_xlabel('horizontal distance [km]')
            ax.set_ylabel('height AGL [m]')
            # ax.set_yticks(np.arange(0,len(plume.lvl),10))
            # ax.set_yticklabels(plume.lvl[::10])
            axh.clear()
            axh.set_ylim([0,150])
            axh.set_xlim([0,dimX*plume.dx])
            # axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
            # axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
            ln = axh.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][n,:], 'r-')
            axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')

            return cntrf, ln, cntr,

        #plot all frames
        ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,cntrf,cntr), interval=3)
        # plt.show()
        ani.save(plume.figdir + 'anim/u/u%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        plt.close()
        print('.....saved in: %s' %(plume.figdir + 'anim/u/u%s.mp4' %Case))
        #
        # #create an animation of vertical velocity---------------------------------------------
        # print('.....creating vertical crossection of W + PM2.5 animation')
        # fig = plt.figure(figsize=(8,6))
        # ax = plt.gca()
        # ax.set_title('VERTICAL VELOCITY : %s' %Case)
        # maxW = int(np.max(csdict['w']))
        # maxPM = int(np.max(csdict['pm25']))
        # wLevels = np.arange(-maxW,maxW+.1,maxW/20)
        # pmLevels = np.arange(maxPM*0.05,maxPM,maxPM/40)
        # # create initial frame
        # # ---u contours and colorbar
        # cntrf = ax.contourf(csdict['w'][0,:,:], cmap=plt.cm.PRGn_r,levels=wLevels,extend='both')
        # cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
        # cbarf.set_label('vertical velocity $[m s^{-1}]$')
        # ax.set_xlabel('horizontal distance [km]')
        # ax.set_ylabel('height AGL [m]')
        # ax.set_xlim([0,dimX])
        # # ---non-filled vapor contours and colorbar
        # # cntr = ax.contour(csdict['pm25'][0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
        # cntr = ax.contour(csdict['pm25'][0,:,:], cmap=plt.cm.Greys,linewidths=2)
        # #
        # # ---heat flux
        # axh = ax.twinx()
        # axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        # axh.set_ylim([0,140])
        # axh.set_xlim([0,dimX])
        # axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        # axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        # axh.tick_params(axis='y', colors='red')
        # ln = axh.plot(csdict['ghfx'][0,:], 'r-')
        # # fig.tight_layout()
        #
        #
        # def update_plot(n,csdict,cntrf,cntr):
        #     ax.clear()
        #     ax.set_title('VERTICAL VELOCITY: %s' %Case)
        #     cntrf = ax.contourf(csdict['w'][n,:,:],cmap=plt.cm.PRGn_r, levels=wLevels,extend='both')
        #     cntr = ax.contour(csdict['pm25'][n,:,:], cmap=plt.cm.Greys,levels=pmLevels,linewidths=0.6)
        #     ax.set_xlabel('horizontal distance [km]')
        #     ax.set_ylabel('height AGL [m]')
        #     ax.set_yticks(np.arange(0,len(plume.lvl),10))
        #     ax.set_yticklabels(plume.lvl[::10])
        #     axh.clear()
        #     axh.set_ylim([0,140])
        #     axh.set_xlim([0,dimX])
        #     axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        #     axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        #     ln = axh.plot(csdict['ghfx'][n,:], 'r-')
        #     axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        #     # fig   .tight_layout()
        #     return cntrf, ln, cntr,
        #
        # #plot all frames
        # ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,cntrf,cntr), interval=3)
        # # plt.show()
        # ani.save(plume.figdir + 'anim/w/w%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        # plt.close()
        # print('.....saved in: %s' %(plume.figdir + 'anim/w/w%s.mp4' %Case))

        #
        # #create an animation of temperature anomaly velocity---------------------------------------------
        # print('.....creating vertical crossection of delT + PM2.5 animation')
        # fig = plt.figure(figsize=(8,6))
        # ax = plt.gca()
        # ax.set_title('TEMPERATURE ANOMALY: %s' %Case)
        # maxT = int(np.max(csdict['temp']))
        # maxPM = int(np.max(csdict['pm25']))
        # uLevels = np.arange(-maxT,maxT+.1,maxT/20.)
        # pmLevels = np.arange(maxPM*0.05,maxPM,maxPM/40.)
        # # create initial frame
        # # ---T contours and colorbar
        # cntrf = ax.contourf(csdict['u'][0,:,:], cmap=plt.cm.Spectral_r,levels=uLevels,extend='both')
        # cbarf = fig.colorbar(cntrf, orientation='horizontal',fraction=0.046, pad=0.1)
        # cbarf.set_label('horizontal velocity $[m s^{-1}]$')
        # ax.set_xlabel('horizontal distance [km]')
        # ax.set_ylabel('height AGL [m]')
        # ax.set_xlim([0,dimX])
        # # ---non-filled vapor contours and colorbar
        # # cntr = ax.contour(csdict['pm25'][0,:,:], cmap=plt.cm.Greys,levels=np.arange(0,2.1,0.3),linewidths=2)
        # cntr = ax.contour(csdict['pm25'][0,:,:], cmap=plt.cm.Greys,linewidths=2)
        #
        # #
        # # ---heat flux
        # axh = ax.twinx()
        # axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        # axh.set_ylim([0,140])
        # axh.set_xlim([0,dimX])
        # axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        # axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        # axh.tick_params(axis='y', colors='red')
        # ln = axh.plot(csdict['ghfx'][0,:], 'r-')
        # # fig.tight_layout()
        #
        #
        # def update_plot(n,csdict,cntrf,cntr):
        #     ax.clear()
        #     ax.set_title('HORIZONTAL VELOCITY CONVERGENCE: %s' %Case)
        #     cntrf = ax.contourf(csdict['u'][n,:,:],cmap=plt.cm.Spectral_r, levels=uLevels,extend='both')
        #     cntr = ax.contour(csdict['pm25'][n,:,:], cmap=plt.cm.Greys,levels=pmLevels,linewidths=0.6)
        #     ax.set_xlabel('horizontal distance [km]')
        #     ax.set_ylabel('height AGL [m]')
        #     ax.set_yticks(np.arange(0,len(plume.lvl),10))
        #     ax.set_yticklabels(plume.lvl[::10])
        #     axh.clear()
        #     axh.set_ylim([0,140])
        #     axh.set_xlim([0,dimX])
        #     axh.set_xticks(np.arange(0,dimX,int(dimX/10)))
        #     axh.set_xticklabels((np.arange(0,dimX,int(dimX/10))*plume.dx/1000).astype(int))
        #     ln = axh.plot(csdict['ghfx'][n,:], 'r-')
        #     axh.set_ylabel('ground heat flux $[kW m^{-2}]$', color='r')
        #     # fig   .tight_layout()
        #
        #     return cntrf, ln, cntr,
        #
        # #plot all frames
        # ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,cntrf,cntr), interval=3)
        # # plt.show()
        # ani.save(plume.figdir + 'anim/u/u%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        # plt.close()
        # print('.....saved in: %s' %(plume.figdir + 'anim/u/u%s.mp4' %Case))
