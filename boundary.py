import numpy as np
import loop
import speed

def apply(i,state) :
    zM = (state.z[i+1,:] > state.zmax)
    zm = (state.z[i+1,:] < state.zmin)
    indM = np.where(zM)[0]
    indm = np.where(zm)[0]

    state.bdy_bot[indM] += 1
    state.bdy_top[indm] += 1

    #Bottom boundary
    if indM.size > 0 :

        ds = (state.zmax - state.z[i,zM]) / (state.C[i,zM]*state.Y[i,zM])
        loop.ray_step(i,zM,ds,state)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],1,0, state.s_dim)

        #Ray normal and tangent
        nx = -state.C[i+1,zM]*state.Y[i+1,zM]
        nz = state.C[i+1,zM]*state.X[i+1,zM]
        tx = state.C[i+1,zM]*state.X[i+1,zM]
        tz = state.C[i+1,zM]*state.Y[i+1,zM]

        #z > 0 upward
        alpha = -tz #t_ray * boundary normal
        beta = -tx #t_ray * boundary tangent (right handed coordinate system)

        cn = dcdz * nz + dcdr * nx
        cs = dcdz * tz + dcdr * tx

        M = beta/alpha
        N = M * (4*cn - 2*M*cs)/state.C[i+1,zM]**2

        state.q[i+1,zM] = state.q[i+1,zM]
        state.p[i+1,zM] = state.p[i+1,zM] + state.q[i+1,zM] * N
        state.Y[i+1,zM] = -state.Y[i+1,zM]


    #Top boundary
    if indm.size > 0 :

        ds = (state.zmin - state.z[i,zm]) / (state.C[i,zm]*state.Y[i,zm])
        loop.ray_step(i,zm,ds,state)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zm],state.r[i+1,zm],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zm],state.r[i+1,zm],1,0, state.s_dim)

        #Ray incident normal and tangent
        nx = -state.C[i+1,zm]*state.Y[i+1,zm]
        nz = state.C[i+1,zm]*state.X[i+1,zm]
        tx = state.C[i+1,zm]*state.X[i+1,zm]
        tz = state.C[i+1,zm]*state.Y[i+1,zm]

        alpha = tz #t_ray * boundary normal
        beta = tx #t_ray * boundary tangent

        cn = dcdz * nz + dcdr * nx
        cs = dcdz * tz + dcdr * tx

        M = beta/alpha
        N = M * (4*cn - 2*M*cs)/state.C[i+1,zm]**2

        state.q[i+1,zm] = state.q[i+1,zm]
        state.p[i+1,zm] = state.p[i+1,zm] + state.q[i+1,zm] * N

        state.Y[i+1,zm] = -state.Y[i+1,zm]
