import numpy as np
import loop
import speed

def apply(i,state) :

    if state.rd_bathy == 0 :
        z_bot = state.zmax
        zM = (state.z[i+1,:] > z_bot)

        nz_bt_bdy = -1 * np.ones_like(state.z[i+1,:])
        tx_bt_bdy = -1 * np.ones_like(state.z[i+1,:])
        nx_bt_bdy = np.zeros_like(state.z[i+1,:])
        tz_bt_bdy = np.zeros_like(state.z[i+1,:])

    elif state.rd_bathy == 1 :

        #array (1,nr) avec les points min et max de la bathy
        bthy_m = np.argmin(state.zmax_r[:,None] < state.r[i+1,:],axis = 0)-1
        bthy_M = bthy_m + 1 #(state.zmax_r[:,None] > state.r[i+1,:])

        dzmax = state.zmax[bthy_M] - state.zmax[bthy_m]
        drmax = state.zmax_r[bthy_M] - state.zmax_r[bthy_m]

        alpha_r = -np.arctan(dzmax/drmax)

        #normal and tangent to the bathy section
        nx_bt_bdy = np.sin(alpha_r)
        nz_bt_bdy = np.cos(alpha_r) #z > 0 upward
        tx_bt_bdy = nz_bt_bdy
        tz_bt_bdy = -nx_bt_bdy

        #depth at the point (interpolate between bthy_m and bthy_M)
        A = np.abs(state.zmax_r[bthy_M] - state.r[i+1,:])
        B = np.abs(state.zmax_r[bthy_M] - state.zmax_r[bthy_m])
        C = np.abs(state.zmax[bthy_m] - state.zmax[bthy_M])

        z_bot = np.max(np.array([state.zmax[bthy_M],state.zmax[bthy_m]])) - C*A/B #marche pô

        #wikipedia linear interpolation
        z_bot = (state.zmax[bthy_m] * ( state.zmax_r[bthy_M] - state.r[i+1,:] ) + state.zmax[bthy_M] * ( state.r[i+1,:] - state.zmax_r[bthy_m] ) )/ B

        r_out = state.r[i+1,:] > state.rmax
        z_bot[r_out] = np.max(state.zmax)

        zM = (state.z[i+1,:] > z_bot)

    else :
        raise ValueError('Bad bathy option')

    zm = (state.z[i+1,:] < state.zmin)
    indM = np.where(zM)[0]
    indm = np.where(zm)[0]

    state.bdy_bot[indM] += 1
    state.bdy_top[indm] += 1

    #Bottom boundary
    if indM.size > 0 :

        if state.rd_bathy == 0 :
            ds = (z_bot - state.z[i,zM]) / (state.C[i,zM]*state.Y[i,zM])
        elif state.rd_bathy == 1 :
            ds = recalculate_step(state, i, z_bot, zM, bthy_m, bthy_M, nx_bt_bdy, nz_bt_bdy)


        loop.ray_step(i,zM,ds,state)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],1,0, state.s_dim)

        #Ray normal and tangent
        nx = -state.C[i+1,zM]*state.Y[i+1,zM]
        nz = state.C[i+1,zM]*state.X[i+1,zM]
        tx = state.C[i+1,zM]*state.X[i+1,zM]
        tz = state.C[i+1,zM]*state.Y[i+1,zM]

        #z > 0 upward
        alpha = tx*nx_bt_bdy[zM] + nz_bt_bdy[zM]*tz #t_ray * boundary normal
        beta = tx_bt_bdy[zM]*tx + tz_bt_bdy[zM]*tz #t_ray * boundary tangent (right handed coordinate system)

        cn = dcdz * nz + dcdr * nx
        cs = dcdz * tz + dcdr * tx

        M = beta/alpha
        N = M * (4*cn - 2*M*cs)/state.C[i+1,zM]**2

        state.q[i+1,zM] = state.q[i+1,zM]
        state.p[i+1,zM] = state.p[i+1,zM] + state.q[i+1,zM] * N

        state.Y[i+1,zM] = state.Y[i+1,zM] - 2 * alpha * nz_bt_bdy[zM] / state.C[i+1,zM]
        state.X[i+1,zM] = state.X[i+1,zM] - 2 * alpha * nx_bt_bdy[zM] / state.C[i+1,zM]

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


def recalculate_step(state, i, z_bot, zM, bthy_m, bthy_M, nx_bt_bdy, nz_bt_bdy) :

    d0_r = state.r[i,:] - state.zmax_r[bthy_m]
    d0_z = state.z[i,:] - state.zmax[bthy_m]

    de_0 = d0_r * nx_bt_bdy + d0_z * nz_bt_bdy

    tx = state.C[i,:]*state.X[i,:]
    tz = state.C[i,:]*state.Y[i,:]

    ds = - de_0 / (tx*nx_bt_bdy + tz*nz_bt_bdy)

    return ds[zM]
