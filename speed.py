import numpy as np
import scipy.stats as stats
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import netCDF4 as nc
#import theo_modes
#import eddy

plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"],
        "font.weight": "bold",
        "text.latex.preamble": [r'\boldmath']})
plt.rc('axes', linewidth=2)

def Munk_c(z,x):
        #Munk profile, from Sensitivity of travel times, Smirnov

        z,x = np.meshgrid(z,x)

        c0 = 1490 #m/s
        za = 1000 #m
        eps = 0.0057
        B = 1000 #m
        L = 1000 #m

        eta = 2*(z - za)/B
        c = c0*(1+eps*(np.exp(-eta) + eta - 1))

        c = c.T
    #print('munk : ', np.shape(c))

        return c#c0 * np.ones_like(z)

#https://stackoverflow.com/questions/20357854/how-can-i-generate-gaussian-random-process-using-matlab
def randg(N,rL,h,cl) :
        x = np.linspace(-rL/2,rL/2,N)
        Z = np.random.normal(loc = 0, scale = h, size=(N))
        F = np.exp(-x**2/(2*cl**2))
        f = np.sqrt(2/np.sqrt(np.pi)) * np.sqrt(rL/N/cl) * np.fft.ifft(np.fft.fft(Z) * np.fft.fft(F))
        return f

def GenerateRandom(z,x, state) :

    z_n = len(z)
    x_n = len(x)
    #print(z_n, x_n)
    #z_,x_ = np.meshgrid(z,x,indexing='ij')

    if state.load_c == 1 :
            #ncf=nc.Dataset('/home6/lheral/Documents/Stage_M2/rec_cel_n3_m.nc')
            ncf=nc.Dataset('../rec_cel_sim_m.nc')

            if state.mean_prof == 0 :
                    c = ncf.variables['C_rec'][...]
                    rn = state.exp_ind
                    c = c[rn,:].T
            else :
                    c = ncf.variables['C_m'][...]

            ncf.close()

    else :
            c = Munk_c(z,x)

    #if state.s_dim == 1 :
            c = c[:,0]



    if state.speed_rand == True :
            #np.random.seed(0)

            #dc = state.speed_std * np.real(randg(z_n, state.zmax, 1, state.speed_Lz)) * np.exp(-z**2/(500**2))

            #dc = theo_modes.solve(z,c)

            dc = eddy.add(x,z)

            c = c + dc

            plt.figure(figsize=(10,8))
            plt.imshow(c, extent=[0,150,0,5000],aspect='auto',cmap='jet',origin='lower')
            cb = plt.colorbar()
            cb.ax.tick_params(labelsize=17)
            cb.set_label('\\textbf{sound speed (m/s)}', rotation=270, fontsize = 17, labelpad = 18)
            plt.ylim(5000,0)
            plt.ylabel('depth (m)')
            plt.xlabel('range (km)')
            plt.xlabel('\\textbf{range (km)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
            plt.ylabel('\\textbf{depth (m)}', fontsize = 17, fontweight = 'bold', labelpad = 10)
            plt.xticks(fontsize = 17)
            plt.yticks(fontsize = 17)
            plt.tick_params(width = 2, length = 4)
            plt.show()


    """

    """

    return c

def Init(z,r,c,z_int,s_dim) :

	#print('shapec',np.shape(c))
	#print(np.shape(z))

        if s_dim == 1 :
                t,cp,k = interp.splrep(z,c,k = 5)
                f_i = interp.BSpline(t,cp,k,extrapolate=True)
                c_int = f_i(z_int)

        elif s_dim == 2 :
                #zr,rr = np.meshgrid(z,r)
                f_i = interp.RectBivariateSpline(z,r,c, kx = 5, ky = 5, s = 0)

        return f_i

def get_speed(z0,r0,f_i, s_dim) :
        if s_dim == 1 :
                z0 = np.array(z0)
                if np.size(z0) > 1:

                        c_int = f_i(z0)

                else :

                        c_int = f_i(z0)

        elif s_dim == 2 :
                z0 = np.array(z0)
                if np.size(z0) > 1:

                        c_int = f_i(z0,r0, grid = False)

                else :

                        c_int = f_i(z0,r0, grid = False)

        #print(c_int)
        return c_int

def get_der(f_i,z0,r0,dz,dr, s_dim) :
        if s_dim == 1 :
                if dr == 0 and dz > 0:
                        df = f_i.derivative(dz)
                        df = df(z0)
                else :
                        df = np.zeros_like(z0)
        elif s_dim == 2 :
                df = f_i.ev(z0,r0, dz, dr)

        return df
