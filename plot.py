import matplotlib.pyplot as plt
import numpy as np

import os
import sys
import numpy as np
from scipy.io import *
import matplotlib.pyplot as plt
sys.path.append ("/home/alexandre/Documents/Stage_M2/Orlando_Python/at/Python/")
#sys.path.append("/home6/lheral/Documents/Stage_M2/Bellhop/Orlando_Python/at/Python/")
from readshd import *
from plotray import *
from read_arrivals_asc import *

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
	"font.weight": "bold",
	"text.latex.preamble": [r'\boldmath']})
plt.rc('axes', linewidth=2)

def show(state) :

    print(np.sum(state.rays_int))
    atten = (state.amp > 0.01)

    if state.exp == "R" :


        #os.system("/home6/lheral/Documents/Stage_M2/Bellhop/at/Bellhop/bellhop.exe ../Bellhop_exp/MunkB_OneBeam")
        os.system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Bellhop/MunkB_Rays")
        fid = open('../Bellhop/MunkB_Rays.bty','r')
        theline = fid.readline()
        theline = fid.readline()
        n       = int( theline )
        rbtykm  = zeros( n )
        zbty    = zeros( n )
        for i in range(n):
            theline = str( fid.readline() )
            datai = theline.split()
            rbtykm[ i] = float( datai[0] )
            zbty[   i] = float( datai[1] )
        fid.close()

        rbty = rbtykm

        #plotray('../Bellhop/MunkB_Rays.ray')
        #plot(rbty,zbty,'m',linewidth=2)


        #plt.figure(figsize = (10,8))
        #system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Code/MunkB_ray")
        #plotray('../Code/MunkB_ray.ray')
        plt.ylim(np.max(state.zmax),0)
        plt.xlim(0,state.rmax/1000)
        #for i in range(state.nr) :
        #    if state.rays_int[i] == False :
        #        continue
        plt.plot(state.r/1000,state.z,'r', linewidth = '1')
        #    print(i)
        plt.plot(state.r/1000,state.z,'ro', linewidth = '1')
        plt.xlabel('\\textbf{Range (km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.ylim(np.max(state.zmax),0)
        plt.xlim(0,state.rmax/1000)

        if state.rd_bathy == 1 :

            plt.plot(state.zmax_r/1000, state.zmax,'b--',linewidth = '2')

            plt.plot(state.zmax_r/1000, state.zmax, 'bo')

            #plt.plot(state.r_c/1000, state.z_c, 'mo')

            #plt.plot(np.array([state.z_c, state.z_c+state.nx_bt_bdy*100])/1000, np.array([state.z_c, state.z_c+state.nz_bt_bdy*100]), 'm-')
            #plt.plot(np.array([state.zmax_r, state.zmax_r+state.tx_bt_bdy[:-1]*100])/1000, np.array([state.zmax, state.zmax+state.tz_bt_bdy[:-1]*100]), 'm-')

            #plt.plot(np.array([state.zmax_r, state.zmax_r+state.nx_node*100])/1000, np.array([state.zmax, state.zmax+state.nz_node*100]), 'k-')
            #plt.plot(np.array([state.zmax_r, state.zmax_r+state.tx_node*100])/1000, np.array([state.zmax, state.zmax+state.tz_node*100]), 'k-')
            #plt.plot(state.r_bot/1000, state.z_bot, 'ko')
            #for i in range(state.nr) :
            #    plt.plot(np.array([state.r[:,i], state.r[:,i]+state.ray_x_bdy[:,i]*100])/1000, np.array([state.z[:,i], state.z[:,i]+state.ray_z_bdy[:,i]*100]), 'm-')
        """
        #plt.grid()

        plt.figure()
        for i in range(state.nr) :
                plt.plot(state.r[:,i],state.ray_x_bdy[:,i], 'k.')
        plt.plot(state.zmax_r, state.nx_node,'r*')
        plt.ylabel('normal x-component')
        plt.xlabel('range')
        """
        plt.figure()
        plt.plot(state.tz)
        """
        plt.figure()
        for i in range(state.nr) :
                plt.plot(state.r[:,i],np.rad2deg(np.arctan(-state.ray_x_bdy[:,i]/state.ray_z_bdy[:,i])), 'k.')
        plt.plot(state.zmax_r, np.rad2deg(np.arctan(-state.nx_node/state.nz_node)),'r*')
        plt.ylabel('normal angle')
        plt.xlabel('range')
        """
        """
        plt.figure()
        plt.plot(state.r/1e3, state.z, 'k')
        plt.ylim(state.zmax,state.zmin)
        plt.xlim(state.rmin,state.rmax/1e3)
        plt.xlabel('range (km)')
        plt.ylabel('depth (m)')
        """

    if state.exp == "TL" :
        """
        system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Code/MunkB_Coh_SGB")

        filename = '../Code/MunkB_Coh_SGB.shd'

        xs = np.nan
        ys = np.nan
        pressure,geometry = readshd(filename,xs,ys)

        zs     = geometry["zs"]
        rarray = geometry["rarray"]; rarraykm = rarray/1000
        zarray = geometry["zarray"]

        Dmax = zarray[-1]
        rmax = rarray[-1]; rmaxkm = rmax/1000

        p_b = np.squeeze( pressure, axis=(0,1) )
        tl = -20*np.log10( abs( p_b ) )

        tl[tl > 120.0] = 120.0

        plt.figure(figsize = (10,8))
        plt.imshow(tl,extent=[0,rmaxkm,0,Dmax],aspect='auto',cmap='jet_r',origin='lower',vmin=60,vmax=120)
        cb = plt.colorbar()
        cb.ax.invert_yaxis()
        cb.ax.tick_params(labelsize=17)
        cb.set_label('\\textbf{TL (dB)}', rotation=270, fontsize = 17, labelpad = 14)

        plt.xlabel('Range (km)')
        plt.ylabel('Depth (m)')
        plt.title('\\textbf{Bellhop}', fontsize = 17, pad = 6)
        plt.xlabel('\\textbf{Range (Km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.ylim(Dmax,0)
        """

        os.system("//home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Bellhop/MunkB_OneBeam")
        filename = '../Bellhop/MunkB_OneBeam.shd'
        fid = open('../Bellhop/MunkB_OneBeam.bty','r')
        theline = fid.readline()
        theline = fid.readline()
        n       = int( theline )
        rbtykm  = zeros( n )
        zbty    = zeros( n )
        for i in range(n):
            theline = str( fid.readline() )
            datai = theline.split()
            rbtykm[ i] = float( datai[0] )
            zbty[   i] = float( datai[1] )
        fid.close()

        rbty = rbtykm

        xs = np.nan
        ys = np.nan
        pressure,geometry = readshd(filename,xs,ys)

        zs     = geometry["zs"]
        rarray = geometry["rarray"]; rarraykm = rarray/1000
        zarray = geometry["zarray"]

        Dmax = zarray[-1]
        rmax = rarray[-1]; rmaxkm = rmax/1000

        p_b = np.squeeze( pressure, axis=(0,1) )
        tl = -20*np.log10( abs( p_b ) )

        tl[tl > 120.0] = 120.0

        print(rbty)


        plt.figure(figsize = (10,8))
        plt.imshow(tl,extent=[0,rmaxkm,0,Dmax],aspect='auto',cmap='jet_r',origin='lower',vmin=40,vmax=120)
        plt.plot(rbty,zbty,'k',linewidth=2)
        cb = plt.colorbar()
        cb.ax.invert_yaxis()
        cb.ax.tick_params(labelsize=17)
        cb.set_label('\\textbf{TL (dB)}', rotation=270, fontsize = 17, labelpad = 14)
        #plt.ylim(np.max(state.zmax),state.zmin)
        #plt.xlim(state.rmin,state.rmax/1000)
        plt.ylim(np.max(state.zmax),0)
        plt.xlim(0,state.rmax/1000)

        zz = np.linspace(state.zmin,np.max(state.zmax),state.Lz)
        rr = np.linspace(state.rmin,state.rmax,state.Lr)
        R, Z = np.meshgrid(rr,zz)

        plt.figure(figsize = (10,8))
        plt.pcolor(R/1000,Z,state.TL, cmap = 'jet_r', vmin = 40, vmax = 120)#, shading = 'nearest')
        #plt.plot(state.r[:,7]/1000, state.z[:,7],'k')
        cbar=plt.colorbar()
        cbar.ax.invert_yaxis()
        cbar.ax.tick_params(labelsize=17)
        cbar.set_label('\\textbf{TL (dB)}', rotation=270, fontsize = 17, labelpad = 14)
        plt.title('\\textbf{FLUOR}', fontsize = 17, pad = 6)
        plt.xlabel('\\textbf{Range (Km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.ylim(np.max(state.zmax),state.zmin)
        plt.xlim(state.rmin,state.rmax/1000)

        if state.rd_bathy == 1 :
            plt.plot(state.zmax_r/1000, state.zmax,'k',linewidth = '4')
            #plt.plot(state.zmax_r/1000, state.zmax, 'ko')

        plt.figure()
        plt.plot(state.r[:,7], state.q[:,7],'k')
        plt.figure()
        plt.plot(state.r[:,7], state.p[:,7],'r')


        plt.figure(figsize = (10,8))
        plt.pcolor(R/1000,Z,state.TL - tl, cmap = 'bwr', vmin = -40, vmax = 40)
        #plt.plot(state.r[:,7]/1000, state.z[:,7],'k')
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=17)
        cbar.set_label('\\textbf{TL difference (dB)}', rotation=270, fontsize = 17, labelpad = 14)
        #plt.title('TL difference')
        plt.ylim(np.max(state.zmax),state.zmin)
        plt.xlim(state.rmin,state.rmax/1000)
        plt.xlabel('\\textbf{Range (Km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        if state.rd_bathy == 1 :
            plt.plot(state.zmax_r/1000, state.zmax,'k',linewidth = '4')


    if state.exp == 'A' :
        """
        os.system("/home6/lheral/Documents/Stage_M2/Bellhop/at/Bellhop/bellhop.exe ../Bellhop_exp/MunkB_Arr")

        fid = open('../Bellhop_exp/MunkB_Arr.bty','r')
        theline = fid.readline()
        theline = fid.readline()
        n       = int( theline )
        rbtykm  = zeros( n )
        zbty    = zeros( n )
        for i in range(n):
            theline = str( fid.readline() )
            datai = theline.split()
            rbtykm[ i] = float( datai[0] )
            zbty[   i] = float( datai[1] )
        fid.close()

        rbty = rbtykm

        Arr, Pos = read_arrivals_asc('../Bellhop/MunkB_Arr.arr', 100 )

        P_ref = 1e-6

        Narr   = np.int_( np.squeeze( Arr['Narr'] ) )
        delay  = np.real( np.squeeze( Arr['delay'] ) )
        A      = abs( np.squeeze( Arr['A'] ) )
        rarray = np.squeeze( Pos['r_range'] )

        plt.figure(figsize = (10,8))
        plt.plot(np.squeeze(Arr['RcvrAngle'])[0:Narr],20*np.log10(A[0:Narr]/P_ref), 'o', label = 'Bellhop')
        plt.plot(state.Angle_rcvr[:state.eigen_ray], 20*np.log10(state.Amp_rcvr[:state.eigen_ray]/P_ref),'r*', label = 'FLUOR')
        plt.xlabel('\\textbf{Angle (°)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Sound pressure level (dB)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.legend(prop={'size': 20, 'weight':'bold'}, loc = 'upper right')
        plt.grid()

        plt.figure(figsize = (10,8))
        plt.plot(np.squeeze(Arr['RcvrAngle'])[0:Narr],np.squeeze(Arr['delay'])[0:Narr], 'o', label = 'Bellhop')
        plt.plot(state.Angle_rcvr[0:state.eigen_ray], state.Delay_rcvr[0:state.eigen_ray],'r*', label = 'FLUOR')
        plt.xlabel('\\textbf{Angle (°)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Delay (s)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.legend(prop={'size': 20, 'weight':'bold'}, loc = 'upper center')
        plt.grid()
        """
        print('eigenrays : ', state.eigen_ray)
        #print('from bellhop : ',Narr)
        """
        os.system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Bellhop/MunkB_OneBeam")
        fid = open('../Bellhop/MunkB_OneBeam.bty','r')
        theline = fid.readline()
        theline = fid.readline()
        n       = int( theline )
        rbtykm  = zeros( n )
        zbty    = zeros( n )
        for i in range(n):
            theline = str( fid.readline() )
            datai = theline.split()
            rbtykm[ i] = float( datai[0] )
            zbty[   i] = float( datai[1] )
        fid.close()

        rbty = rbtykm

        plotray('../Bellhop/MunkB_OneBeam.ray')
        plt.plot(rbty,zbty,'m',linewidth=2)
        plt.ylim(np.max(state.zmax),state.zmin)
        plt.xlim(state.rmin,state.rmax/1000)
        plt.plot(state.r_rcvr/1000,state.z_rcvr,'ro')
        """

        plt.figure()
        for i in range(state.eigen_ray) :
            #plt.plot(state.r[atten[:,i],i]/1000,state.z[atten[:,i],i],'k', linewidth = '1')

            plt.plot(state.r[atten[:,state.ray_num[i]],state.ray_num[i]],state.z[atten[:,state.ray_num[i]],state.ray_num[i]],'k')
        plt.ylim(np.max(state.zmax),state.zmin)
        plt.xlim(state.rmin,state.rmax)
        plt.plot(state.zmax_r, state.zmax,'r',linewidth = '2')
        plt.plot(state.r_rcvr,state.z_rcvr,'ro')

        print(state.eigen_ray)

        plt.figure()
        plt.plot(state.r[:,state.ray_num],state.z[:,state.ray_num],'k', linewidth = '1')
        plt.plot(state.r_rcvr,state.z_rcvr,'ro')
        plt.ylim(np.max(state.zmax),state.zmin)
        plt.xlim(state.rmin,state.rmax)
        plt.plot(state.zmax_r, state.zmax,'r',linewidth = '2')


        plt.figure()
        plt.plot(state.Angle_rcvr[:state.eigen_ray], 20*np.log10(state.Amp_rcvr[:state.eigen_ray]/1e-6),'r*')
        plt.title("amp db")
        plt.xlabel('rcvr angle')

        plt.figure()
        plt.plot(state.Angle_rcvr[0:state.eigen_ray], state.Delay_rcvr[0:state.eigen_ray],'r*')
        plt.title("delay")
        plt.xlabel('rcvr angle')


    plt.show()
