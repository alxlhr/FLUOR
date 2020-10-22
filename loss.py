from numba import jit
import numpy as np
#import fortran_loss as floss
import matplotlib.pyplot as plt

def calc_normals(state) :
    state.tx = state.C*state.X
    state.tz = state.C*state.Y
    state.nx = -state.C*state.Y
    state.nz = state.C*state.X


def calc_TL(state) :
    #state.dangle0 = np.zeros_like(state.q)
    #state.dangle0[:,1:] = (state.angle[:,1:] - state.angle[:,:-1])
    state.W = np.abs(state.q * state.dangle0[None,:])
    state.P = np.zeros((state.Lz,state.Lr))
    zz = np.linspace(state.zmin,np.max(state.zmax),state.Lz)
    rr = np.linspace(state.rmin,state.rmax,state.Lr)
    R, Z = np.meshgrid(rr,zz)

    print('shape P', np.shape(state.P))

    if state.use_fortran == 1 :
        #floss.loop(state.W, R, Z, state.r, state.z, state.nx, state.nz, state.tz, state.T, state.q, state.P)
        state.P = floss.loop(q = state.q, t = state.T, tz = state.tz, nz = state.nz, nx = state.nx, z = state.z, r = state.r, zz = Z, rr = R, w = state.W, lr = state.Lr, lz = state.Lz, nr = state.nr, nmax = state.n_max, f = state.f)
    else :
        loop(state)

    state.TL = -20*np.log10(4*np.pi*np.abs(state.P))
    state.TL[np.isnan(state.TL)] = 120
    state.TL[np.isinf(state.TL)] = 120


#@jit(nopython = True)#, parallel = True)
def loop(state) :
    zz = np.linspace(state.zmin,np.max(state.zmax),state.Lz)
    rr = np.linspace(state.rmin,state.rmax,state.Lr)
    R, Z = np.meshgrid(rr,zz)
    rec = np.zeros_like(R, dtype = bool)
    n_ = np.zeros_like(R)
    #s_ = np.zeros_like(R)
    T_ = np.zeros_like(R)
    q_ = np.zeros_like(R)
    al = np.zeros_like(R)
    #W_ = np.zeros_like(R)
    r_ = np.zeros_like(R)
    n_val = np.zeros_like(R)
    A_ = np.zeros_like(R, dtype=complex)
    print(state.rays_int)
    j = 7
    #for j in range(state.nr) :
    print("beam %i / %i" %(j, state.nr), end = '\r')
    if (state.rays_int[j] == True) :
        for i in range(2,state.n_max-1) :
            if state.W[i,j] > 0:
                #coords influenced by the ray j between i and i+1
                rec[:,:] = ((R >= state.r[i,j]) & (R < state.r[i+1,j]))
                #normal distance all the points and segment i

                n_[:,:] = np.abs( (R - state.r[i,j])*state.nx[i,j] + (Z - state.z[i,j])*state.nz[i,j])

                #print(np.shape(n_))
                #along ray coordinate
                #s_[:,:] = (R - state.r[i-1,j])*state.tx[i-1,j] + (Z - state.z[i-1,j])*state.tz[i-1,j]

                al[:,:] = ((R - state.r[i-1,j])*state.tx[i-1,j] + (Z - state.z[i-1,j])*state.tz[i-1,j]) / np.sqrt((state.r[i-1,j] - state.r[i,j])**2 + (state.z[i-1,j] - state.z[i,j])**2)

                T_[:,:] = state.T[i-1,j] + al * (state.T[i,j] - state.T[i-1,j])
                q_[:,:] = state.q[i-1,j] + al * (state.q[i,j] - state.q[i-1,j])
                W_ = state.W[i,j] #+ al * (state.W[i,j] - state.W[i-1,j])
                r_[:,:] = state.r[i-1,j] + al * (state.r[i,j] - state.r[i-1,j])

                n_val[:,:] = (W_ - n_*rec)/W_*rec
                n_val[n_val < 0] = 0

                A_[:,:] = state.amp[i,j] * 1/(4*np.pi) * (-1j)**state.m[i,j] * np.sqrt(np.abs(state.C[i,j]*np.cos(state.angle_0[j])/(r_*state.c0*q_)))
                """
                print(np.shape(n_[rec]))
                "
                if i > 50 :
                    plt.figure()
                    plt.pcolormesh(R,Z, n_val * A_ * np.exp(1j*2*np.pi*state.f*T_)*np.exp(1j*state.phi[i,j]), cmap = 'jet')
                    plt.colorbar()
                    plt.plot(state.r[:,j],state.z[:,j])
                    plt.plot(state.r[:,j],state.z[:,j],'k.')
                    plt.plot(state.r[i,j], state.z[i,j],'ro')
                    plt.ylim(200,0)
                    plt.show()
                """
                
                state.P = state.P + n_val * A_ * np.exp(1j*2*np.pi*state.f*T_)*np.exp(1j*state.phi[i,j])
        """

def loop(state) :
    zz = np.linspace(state.zmin,np.max(state.zmax),state.Lz)
    rr = np.linspace(state.rmin,state.rmax,state.Lr)
    R, Z = np.meshgrid(rr,zz)
    #j = 7
    for j in range(state.nr) :
        print("beam %i / %i" %(j, state.nr), end = '\r')

        rec = ((R[:,:,None] >= state.r[:-1,j]) & (R[:,:,None] < state.r[1:,j]))
        #Receivers in ray-centered coordinates
        n_ = np.abs( (R[:,:,None] - state.r[1:,j])*state.nx[1:,j] + (Z[:,:,None] - state.z[1:,j])*state.nz[1:,j])
        s_ = (R[:,:,None] - state.r[:-1,j])*state.tx[:-1,j] + (Z[:,:,None] - state.z[:-1,j])*state.tz[:-1,j]

        #Proportional distance along the ray
        al = s_ / np.sqrt((state.r[:-1,j] - state.r[1:,j])**2 + (state.z[:-1,j] - state.z[1:,j])**2)

        T_ = state.T[:-1,j] + al * (state.T[1:,j] - state.T[:-1,j])
        q_ = state.q[:-1,j] + al * (state.q[1:,j] - state.q[:-1,j])
        W_ = state.W[:-1,j] + al * (state.W[1:,j] - state.W[:-1,j])
        r_ = state.r[:-1,j] + al * (state.r[1:,j] - state.r[:-1,j])

        n_val = (n_ <= W_)

        plt.figure()
        plt.pcolormesh(n_[:,:,0] <= W_[:,:,0], cmap = 'jet')
        plt.colorbar()
        plt.show()


        A_ = state.amp[1:,j] * 1/(4*np.pi) * (-1j)**state.m[1:,j] * np.sqrt(np.abs(state.C[1:,j]*np.cos(state.angle_0[j])/(r_*state.c0*q_)))

        state.P = state.P + np.sum((W_ - n_)/W_*n_val*rec * A_ * np.exp(1j*2*np.pi*state.f*T_)*np.exp(1j*state.phi[1:,j]), axis = 2)

        print(np.shape(state.P))


        for i in range(1,state.n_max) :
            if state.W[i,j] > 0 :
                #coords influenced by the ray j between i and i+1
                rec = ((R >= state.r[i-1,j]) & (R < state.r[i,j]))
                #normal distance between those points and segment i
                n_ = np.abs( (R - state.r[i,j])*state.nx[i,j] + (Z - state.z[i,j])*state.nz[i,j])
                #along ray coordinate
                s_ = (R - state.r[i-1,j])*state.tx[i-1,j] + (Z - state.z[i-1,j])*state.tz[i-1,j]

                al = s_ / np.sqrt((state.r[i-1,j] - state.r[i,j])**2 + (state.z[i-1,j] - state.z[i,j])**2)

                T_ = state.T[i-1,j] + al * (state.T[i,j] - state.T[i-1,j])
                q_ = state.q[i-1,j] + al * (state.q[i,j] - state.q[i-1,j])
                W_ = state.W[i-1,j] + al * (state.W[i,j] - state.W[i-1,j])
                r_ = state.r[i-1,j] + al * (state.r[i,j] - state.r[i-1,j])

                n_val = (n_ <= W_)

                A_ = state.amp[i,j] * 1/(4*np.pi) * (-1j)**state.m[i,j] * np.sqrt(np.abs(state.C[i,j]*np.cos(state.angle_0[j])/(r_*state.c0*q_)))

                state.P = state.P + (W_ - n_)/W_*n_val*rec * A_ * np.exp(1j*2*np.pi*state.f*T_)*np.exp(1j*state.phi[i,j])

"""
