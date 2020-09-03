from numba import jit
import numpy as np

def calc_normals(state) :
    state.tx = state.C*state.X
    state.tz = state.C*state.Y
    state.nx = -state.C*state.Y
    state.nz = state.C*state.X


def calc_TL(state) :
    state.W = np.abs(state.q * state.dangle0[np.newaxis,:])
    state.P = np.zeros((state.Lz,state.Lr))

    loop(state)

    state.TL = -20*np.log10(4*np.pi*np.abs(state.P))
    state.TL[np.isnan(state.TL)] = 120
    state.TL[np.isinf(state.TL)] = 120


#@jit(nopython = True)#, parallel = True)
def loop(state) :
    zz = np.linspace(state.zmin,state.zmax,state.Lz)
    rr = np.linspace(state.rmin,state.rmax,state.Lr)
    R, Z = np.meshgrid(rr,zz)
    #j = 7
    for j in range(state.nr) :
        print("beam %i / %i" %(j, state.nr), end = '\r')
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

                A_ = 1/(4*np.pi) * (-1j)**state.m[i,j] * np.sqrt(np.abs(state.C[i,j]*np.cos(state.angle_0[j])/(r_*state.c0*q_)))

                state.P = state.P + (W_ - n_)/W_*n_val*rec * A_ * np.exp(1j*2*np.pi*state.f*T_)
