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
    csdict = plume.load_CS_prep_Profiles(Case)
    T0 = np.load(plume.wrfdir + 'interp/profT0' + Case + '.npy')
    U0 = np.load(plume.wrfdir + 'interp/profU0' + Case + '.npy')

    dimT, dimZ, dimX = np.shape(csdict['temp'])
    dimY = np.shape(csdict['ghfx2D'])[1]
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

    if doAnimations:
        maxPM = int(np.max(csdict['pm25']))
        pmLevels = np.geomspace(30,maxPM/10,num=10)
        PMcontours = ma.masked_where(csdict['pm25'] <= 30,csdict['pm25'])

        # #create an animation of horizontal velocity---------------------------------------------
        # fig = plt.figure(figsize=(22,4))
        # gs = fig.add_gridspec(ncols=2, nrows=1,width_ratios=[6,1])
        # print('.....creating vertical crossection of U + PM2.5 animation')
        # plt.suptitle('RELATIVE HORIZONTAL VELOCITY: %s' %Case)
        #
        # ax1=fig.add_subplot(gs[0])
        # axh1=ax1.twinx()
        #
        # # create initial frame
        # # ---u contours and colorbar
        # im = ax1.imshow((csdict['u'][0,:,:].T - U0).T, origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)
        # cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
        # cbari.set_label('ralative horizontal velocity $[m s^{-1}]$')
        # # ---non-filled vapor contours and colorbar
        # cntr = ax1.contour(PMcontours[0,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
        # ax1.set_xlabel('horizontal distance [m]')
        # ax1.set_ylabel('height AGL [m]')
        # ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
        #
        # # ---heat flux
        # ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][0,:], 'r-')
        # axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
        # axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])
        # axh1.tick_params(axis='y', colors='red')
        #
        #
        # ax2=fig.add_subplot(gs[1])
        # fim = ax2.imshow(csdict['ghfx2D'][0,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
        # ax2.set_aspect('equal')
        # cbarif = fig.colorbar(fim, orientation='horizontal')
        # cbarif.set_label('heat flux [$kW / m^2$]')
        # ax2.set_xlabel('x distance [m]')
        # ax2.set_ylabel('y distance [m]')
        # plt.tight_layout(rect=[0, 0, 1, 0.95])
        #
        # def update_plot(n,csdict,cntrf,cntr):
        #     ax1.clear()
        #     im = ax1.imshow((csdict['u'][n,:,:].T - U0).T,origin='lower',extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.RdBu_r,vmin=-4, vmax=4)
        #     cntr = ax1.contour(PMcontours[n,:,:], extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.7)
        #     ax1.set_xlabel('horizontal distance [m]')
        #     ax1.set_ylabel('height AGL [m]')
        #     ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
        #
        #     axh1.clear()
        #     axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])
        #     ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][n,:], 'r-')
        #     axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
        #
        #     ax2.clear()
        #     fim = ax2.imshow(csdict['ghfx2D'][n,75:175,0:75],cmap=plt.cm.YlOrRd,extent=[0,3000,3000,7000],vmin=0, vmax = 150)
        #     ax2.set_aspect('equal')
        #     ax2.set_xlabel('x distance [m]')
        #     ax2.set_ylabel('y distance [m]')
        #
        #     return cntrf, ln, cntr,
        #
        # #plot all frames
        # ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,im,cntr), interval=3)
        # # plt.show()
        # ani.save(plume.figdir + 'anim/u/u%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        # plt.close()
        # print('.....saved in: %s' %(plume.figdir + 'anim/u/u%s.mp4' %Case))


        #create an animation of vertical velocity---------------------------------------------
        fig = plt.figure(figsize=(22,4))
        gs = fig.add_gridspec(ncols=2, nrows=1,width_ratios=[6,1])

        print('.....creating vertical crossection of W + PM2.5 animation')
        plt.suptitle('RELATIVE VERTICAL VELOCITY: %s' %Case)

        ax1=fig.add_subplot(gs[0])
        # create initial frame
        # ---w contours and colorbar
        im = ax1.imshow(csdict['w'][0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
        cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
        cbari.set_label('vertical velocity $[m s^{-1}]$')
        # ---non-filled vapor contours and colorbar
        cntr = ax1.contour(PMcontours[0,:,:],extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.7)
        ax1.set_xlabel('horizontal distance [m]')
        ax1.set_ylabel('height AGL [m]')
        ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')
        # ---heat flux
        axh1 = ax1.twinx()
        ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][0,:], 'r-')
        axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
        axh1.tick_params(axis='y', colors='red')
        axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])

        ax2=fig.add_subplot(gs[1])
        fim = ax2.imshow(csdict['ghfx2D'][0,75:175,0:75],cmap=plt.cm.YlOrRd, extent=[0,3000,3000,7000],vmin=0, vmax = 150)
        ax2.set_aspect('equal')
        cbarif = fig.colorbar(fim, orientation='horizontal')
        cbarif.set_label('heat flux [$kW / m^2$]')
        ax2.set_xlabel('x distance [m]')
        ax2.set_ylabel('y distance [m]')
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        def update_plot(n,csdict,cntrf,cntr):
            ax1.clear()
            im = ax1.imshow((csdict['u'][n,:,:].T - U0).T,origin='lower',extent=[0,dimX*plume.dx,0,plume.lvl[-1]],cmap=plt.cm.PRGn_r, vmin=-5, vmax=5)
            cntr = ax1.contour(PMcontours[n,:,:], extent=[0,dimX*plume.dx,0,plume.lvl[-1]],levels=pmLevels,locator=ticker.LogLocator(),cmap=plt.cm.Greys,linewidths=0.6)
            ax1.set_xlabel('horizontal distance [m]')
            ax1.set_ylabel('height AGL [m]')
            ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal')

            axh1.clear()
            ln = axh1.plot(np.arange(dimX)*plume.dx, csdict['ghfx'][n,:], 'r-')
            axh1.set_ylabel('fire heat flux $[kW m^{-2}]$', color='r')
            axh1.set(xlim=[0,dimX*plume.dx],ylim=[0,150])

            ax2.clear()
            fim = ax2.imshow(csdict['ghfx2D'][n,75:175,0:75],cmap=plt.cm.YlOrRd,extent=[0,3000,3000,7000],vmin=0, vmax = 150)
            ax2.set_aspect('equal')
            ax2.set_xlabel('x distance [m]')
            ax2.set_ylabel('y distance [m]')

            return cntrf, ln, cntr,

        #plot all frames
        ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,im,cntr), interval=3)
        # plt.show()
        ani.save(plume.figdir + 'anim/w/w%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
        plt.close()
        print('.....saved in: %s' %(plume.figdir + 'anim/w/w%s.mp4' %Case))

        if Case == 'W4F7R4':
            fig = plt.figure(figsize=(12,4))

            ax1 = plt.gca()
            # create initial frame
            im = ax1.imshow(csdict['pm25'][0,:,:], origin='lower', extent=[0,dimX*plume.dx,0,plume.lvl[-1]],vmin = 0, vmax=np.max(csdict['pm25'][-1,:,:])/12.,cmap=plt.cm.cubehelix_r)
            cbari = fig.colorbar(im, orientation='horizontal',aspect=60, shrink=0.5)
            cbari.set_label('concentration [ug/kg]')
            ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal',xlabel='horizontal distance [m]', ylabel='height AGL [m]')
            plt.tight_layout()
            def update_plot(n,csdict,cntrf,cntr):
                ax1.clear()
                im = ax1.imshow(csdict['pm25'][n,:,:],origin='lower',extent=[0,dimX*plume.dx,0,plume.lvl[-1]],vmin = 0, vmax=np.max(csdict['pm25'][-1,:,:])/12.,cmap=plt.cm.cubehelix_r)
                ax1.set(xlim=[0,dimX*plume.dx],ylim=[0,plume.lvl[-1]],aspect='equal', xlabel='horizontal distance [m]', ylabel='height AGL [m]')
                return cntrf, ln, cntr,

            #plot all frames
            ani=animation.FuncAnimation(fig, update_plot, dimT, fargs=(csdict,im,cntr), interval=3)
            # plt.show()
            ani.save(plume.figdir + 'anim/w/smoke%s.mp4' %Case, writer='ffmpeg',fps=10, dpi=250)
            plt.close()

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
