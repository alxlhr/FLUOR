import numpy as np
import loop
import speed

def apply(i,state) :

    if state.rd_bathy == 0 :
        z_bot = state.zmax
        zM = (state.z[i+1,:] > z_bot)

        nz_bt_bdy = -1
        tx_bt_bdy = -1
        nx_bt_bdy = 0
        tz_bt_bdy = 0

    elif state.rd_bathy == 1 :

        #array (1,nr) avec les points min et max de la bathy
        bthy_m = np.argmin(state.zmax_r[:,None] < state.r[i+1,:],axis = 0)-1
        bthy_M = bthy_m + 1 #(state.zmax_r[:,None] > state.r[i+1,:])

        print(bthy_m)

        #abs for now (la pente est dans la direction de propa du rayon)
        dzmax = np.abs(state.zmax[bthy_M] - state.zmax[bthy_m])
        drmax = np.abs(state.zmax_r[bthy_M] - state.zmax_r[bthy_m])

        alpha_r = np.arctan(dzmax/drmax)

        #normal and tangent to the bathy section
        nx_bt_bdy = np.sin(alpha_r)
        nz_bt_bdy = np.cos(alpha_r)
        tx_bt_bdy = nz_bt_bdy
        tz_bt_bdy = -nx_bt_bdy

        #bottom value
        A = np.abs(state.zmax_r[bthy_M] - state.r[i+1,:])
        B = np.abs(state.zmax_r[bthy_M] - state.zmax_r[bthy_m])
        C = np.abs(state.zmax[bthy_m] - state.zmax[bthy_M])

        z_bot = state.zmax[bthy_M] + C*A/B

        zM = (state.z[i+1,:] > z_bot)

        nz_bt_bdy = -1
        tx_bt_bdy = -1
        nx_bt_bdy = 0
        tz_bt_bdy = 0

    else :
        raise ValueError('Bad bathy option')

    zm = (state.z[i+1,:] < state.zmin)
    indM = np.where(zM)[0]
    indm = np.where(zm)[0]

    state.bdy_bot[indM] += 1
    state.bdy_top[indm] += 1

    #Bottom boundary
    if indM.size > 0 :

        ds = (z_bot[zM] - state.z[i,zM]) / (state.C[i,zM]*state.Y[i,zM])
        loop.ray_step(i,zM,ds,state)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],1,0, state.s_dim)

        #Ray normal and tangent
        nx = -state.C[i+1,zM]*state.Y[i+1,zM]
        nz = state.C[i+1,zM]*state.X[i+1,zM]
        tx = state.C[i+1,zM]*state.X[i+1,zM]
        tz = state.C[i+1,zM]*state.Y[i+1,zM]

        #z > 0 upward
        alpha = tx*nx_bt_bdy + nz_bt_bdy*tz #t_ray * boundary normal
        beta = tx_bt_bdy*tx + tz_bt_bdy*tz #t_ray * boundary tangent (right handed coordinate system)

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
