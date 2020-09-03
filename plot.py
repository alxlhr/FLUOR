import matplotlib.pyplot as plt
import numpy as np

from os import *
import sys
sys.path.append ("/home/alexandre/Documents/Stage_M2/Orlando_Python/at/Python/")
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

    if state.exp == "R" :

        plt.figure(figsize = (10,8))
        #system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Code/MunkB_ray")
        #plotray('../Code/MunkB_ray.ray')
        plt.ylim(5000,0)
        plt.xlim(0,150)
        plt.plot(state.r/1000,state.z,'k', linewidth = '2')
        plt.xlabel('\\textbf{Range (km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)
        plt.grid()
        #plt.title("Bellhop (black) vs Alex (red) - Munk profile")

        """
        plt.figure()
        plt.plot(state.r/1e3, state.z, 'k')
        plt.ylim(state.zmax,state.zmin)
        plt.xlim(state.rmin,state.rmax/1e3)
        plt.xlabel('range (km)')
        plt.ylabel('depth (m)')
        """

    if state.exp == "TL" :

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



        zz = np.linspace(state.zmin,state.zmax,state.Lz)
        rr = np.linspace(state.rmin,state.rmax,state.Lr)
        R, Z = np.meshgrid(rr,zz)

        plt.figure(figsize = (10,8))
        plt.pcolor(R/1000,Z,state.TL, cmap = 'jet_r', vmin = 60, vmax = 120)
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
        plt.ylim(state.zmax,state.zmin)
        plt.xlim(state.rmin,state.rmax/1000)

        plt.figure(figsize = (10,8))
        plt.pcolor(R/1000,Z,state.TL - tl, cmap = 'bwr', vmin = -40, vmax = 40)
        cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize=17)
        cbar.set_label('\\textbf{TL difference (dB)}', rotation=270, fontsize = 17, labelpad = 14)
        #plt.title('TL difference')
        plt.ylim(state.zmax,state.zmin)
        plt.xlim(state.rmin,state.rmax/1000)
        plt.xlabel('\\textbf{Range (Km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.ylabel('\\textbf{Depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tick_params(width = 2, length = 4)

    if state.exp == 'A' :

        system("/home/alexandre/Documents/Stage_M2/Bellhop/at/at/Bellhop/bellhop.exe ../Code/MunkB_Arr")

        Arr, Pos = read_arrivals_asc('../Code/MunkB_Arr.arr', 100 )

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
        print('from bellhop : ',Narr)

        plt.figure()
        for i in range(state.eigen_ray) :
            plt.plot(state.r[:,state.ray_num[i]],state.z[:,state.ray_num[i]],'k')

        plt.plot(state.r_rcvr,state.z_rcvr,'ro')
        plt.ylim(state.zmax,state.zmin)
        plt.xlim(state.rmin,state.rmax)


        plt.figure()
        plt.plot(state.Angle_rcvr[:state.eigen_ray], -10*np.log10(state.Amp_rcvr[:state.eigen_ray]),'r*')
        plt.title("amp db")
        plt.xlabel('rcvr angle')

        plt.figure()
        plt.plot(state.Angle_rcvr[0:state.eigen_ray], state.Delay_rcvr[0:state.eigen_ray],'r*')
        plt.title("delay")
        plt.xlabel('rcvr angle')
        """

    plt.show()
