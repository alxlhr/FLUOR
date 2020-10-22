import numpy as np
import loop
import speed
import bdy_loss

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
        bthy_M = bthy_m + 1

        #depth at the point (interpolate between bthy_m and bthy_M)
        B = np.abs(state.zmax_r[bthy_M] - state.zmax_r[bthy_m])

        #linear interpolation of the depth at that point (linterp from wikipedia)
        z_bot = (state.zmax[bthy_m] * ( state.zmax_r[bthy_M] - state.r[i+1,:] ) + state.zmax[bthy_M] * ( state.r[i+1,:] - state.zmax_r[bthy_m] ) )/ B
        #print(crossing_depth(state.r[i,:],state.z[i,:],state.r[i+1,:],state.z[i+1,:], state.zmax_r[bthy_m], state.zmax[bthy_m], state.zmax_r[bthy_M], state.zmax[bthy_M]))
        
        r_out = (state.r[i+1,:] > state.rmax)
        z_bot[r_out] = state.zmax[-1] #bathy outside the domain to avoid problems with reflections
      
    else :
        raise ValueError('Bad bathy option')

    zM = (state.z[i+1,:] > z_bot)
    zm = (state.z[i+1,:] < state.zmin)
    indM = np.where(zM)[0]
    indm = np.where(zm)[0]

    state.bdy_bot[i+1,:] = state.bdy_bot[i,:] + np.ones((state.nr))*zM
    state.bdy_top[i+1,:] = state.bdy_top[i,:] + np.ones((state.nr))*zm

    #Bottom boundary
    if indM.size > 0 :            
        if state.rd_bathy == 0 :
            ds = (z_bot - state.z[i,zM]) / (state.C[i,zM]*state.Y[i,zM])
        elif state.rd_bathy == 1 :
            nx_bt_bdy = state.nx_bt_bdy[bthy_m[zM]]
            nz_bt_bdy = state.nz_bt_bdy[bthy_m[zM]]
            ds = recalculate_step(state, i, zM, bthy_m,nx_bt_bdy, nz_bt_bdy)

        loop.ray_step(i,zM,ds,state)

        #If linterp, recalculate normals
        if state.bathy_linterp == 1 :
            [nx_bt_bdy, nz_bt_bdy] = interpolate_normals(i, state,zM, bthy_m)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],1,0, state.s_dim)

        tx_bt_bdy = -nz_bt_bdy
        tz_bt_bdy = nx_bt_bdy

        #tangent : c * (X,Y)
        #normal :  c * (-Y,X)

        #z > 0 upward
        alpha = state.tx[i+1,zM]*nx_bt_bdy + nz_bt_bdy*state.tz[i+1,zM] #t_ray * boundary normal
        beta = tx_bt_bdy*state.tx[i+1,zM] + tz_bt_bdy*state.tz[i+1,zM] #t_ray * boundary tangent (right handed coordinate system)

        state.Y[i+1,zM] = state.Y[i+1,zM] - 2 * alpha * nz_bt_bdy / state.C[i+1,zM]
        state.X[i+1,zM] = state.X[i+1,zM] - 2 * alpha * nx_bt_bdy / state.C[i+1,zM]
        
        state.tx[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]
        state.tz[i+1,zM] = state.C[i+1,zM]*state.Y[i+1,zM]
        state.nx[i+1,zM] = -state.C[i+1,zM]*state.Y[i+1,zM]
        state.nz[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]

        #Ray normal and tangent
        theta_I = state.tx[i+1,zM]*nx_bt_bdy + state.tz[i+1,zM]*nz_bt_bdy

        cn = dcdz * state.nz[i+1,zM] + dcdr * state.nx[i+1,zM]
        cs = dcdz * state.tz[i+1,zM] + dcdr * state.tx[i+1,zM]

        M = beta/alpha
        N = M * (4*cn - 2*M*cs)/state.C[i+1,zM]**2

        q_prev = state.q[i+1,zM].copy()
        state.q[i+1,zM] = q_prev
        state.p[i+1,zM] = state.p[i+1,zM] + q_prev * N

        #state.X[i+1,zM] = (- alpha * nx_bt_bdy + beta * tx_bt_bdy) / state.C[i+1,zM]
        #state.Y[i+1,zM] = (- alpha * nz_bt_bdy + beta * tz_bt_bdy) / state.C[i+1,zM]

        #boundary losses
        R = bdy_loss.fluid_fluid(state,i,zM,theta_I)
        state.amp[i+1,zM] = state.amp[i,zM] * np.abs(R)
        state.phi[i+1,zM] = state.phi[i,zM] + np.angle(R)

        #Angles check

        #theta_R = tx*nx_bt_bdy + tz*nz_bt_bdy
        #print("**********")
        #print("I : ", theta_I)
        #print("R : ", theta_R)
        #chck_ang = theta_I + theta_R < 1e-12
        #if (np.any(chck_ang == False)) :
        #    print("Problem with incident/reflected angles : ",i)

    #Top boundary
    if indm.size > 0 :

        ds = (state.zmin - state.z[i,:]) / (state.C[i,:]*state.Y[i,:])
        loop.ray_step(i,zm,ds,state)

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zm],state.r[i+1,zm],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zm],state.r[i+1,zm],1,0, state.s_dim)

        #Ray incident normal and tangent

        alpha = state.tz[i+1,zm] #t_ray * boundary normal
        beta = state.tx[i+1,zm] #t_ray * boundary tangent

        state.Y[i+1,zm] = -state.Y[i+1,zm]

        state.tx[i+1,zm] = state.C[i+1,zm]*state.X[i+1,zm]
        state.tz[i+1,zm] = state.C[i+1,zm]*state.Y[i+1,zm]
        state.nx[i+1,zm] = -state.C[i+1,zm]*state.Y[i+1,zm]
        state.nz[i+1,zm] = state.C[i+1,zm]*state.X[i+1,zm]

        cn = dcdz * state.nz[i+1,zm] + dcdr * state.nx[i+1,zm]
        cs = dcdz * state.tz[i+1,zm] + dcdr * state.tx[i+1,zm]

        M = beta/alpha
        N = M * (4*cn - 2*M*cs)/state.C[i+1,zm]**2

        q_prev = state.q[i+1,zm].copy()
        state.q[i+1,zm] = q_prev
        state.p[i+1,zm] = state.p[i+1,zm] + q_prev * N


def recalculate_step(state, i, zM, bthy_m, nx_bt_bdy, nz_bt_bdy) :

    d0_r = state.r[i,zM] - state.zmax_r[bthy_m[zM]]
    d0_z = state.z[i,zM] - state.zmax[bthy_m[zM]]

    d_r = state.r[i+1,zM] - state.zmax_r[bthy_m[zM]]
    d_z = state.z[i+1,zM] - state.zmax[bthy_m[zM]]

    de_0 = d0_r * nx_bt_bdy + d0_z * nz_bt_bdy
    de = d_r * nx_bt_bdy + d_z * nz_bt_bdy

    tx = state.C[i,zM]*state.X[i,zM]
    tz = state.C[i,zM]*state.Y[i,zM]

    dh = - de_0 / (-de_0 + de) * state.ds0[zM]

    ds = np.zeros_like(state.ds0)
    ds[zM] = - de_0 / (tx*nx_bt_bdy + tz*nz_bt_bdy)

    if (np.any(ds < 0)) :
        print("Problem with negative step size : ",i)
        print(zM)
        state.rays_int[ds < 0] = False

    return ds


def calculate_normals(state) :

    dzmax = np.diff(state.zmax)
    drmax = np.diff(state.zmax_r)

    alpha_r = -np.arctan(dzmax/drmax)

    #normal to the bathy section
    state.nx_bt_bdy = -np.sin(alpha_r)
    state.nz_bt_bdy = -np.cos(alpha_r) #z > 0 upward

    state.tx_bt_bdy = -state.nz_bt_bdy
    state.tz_bt_bdy = state.nx_bt_bdy

    #get center of the bathy sections

    dl = np.sqrt(dzmax**2 + drmax**2)
    state.z_c = dl/2 * np.sin(-alpha_r) + state.zmax[:-1]
    state.r_c =  dl/2 * np.cos(-alpha_r) + state.zmax_r[:-1]

    #Get normals at the nodes

    state.nx_node = np.zeros_like(state.zmax,dtype = float)
    state.nz_node = np.zeros_like(state.zmax,dtype = float)
    state.tx_node = np.zeros_like(state.zmax,dtype = float)
    state.tz_node = np.zeros_like(state.zmax,dtype = float)

    state.nx_node[0] = state.nx_bt_bdy[0]
    state.nz_node[0] = state.nz_bt_bdy[0]
    state.nx_node[-1] = state.nx_bt_bdy[-1]
    state.nz_node[-1] = state.nz_bt_bdy[-1]

    for j in range(1,len(state.zmax)-1) :
        #should be ok :
        dist_bwd = np.sqrt( (state.zmax_r[j] - state.r_c[j-1])**2 + (state.zmax[j] - state.z_c[j-1])**2) #distance from the node to the previous normal
        dist_fwd = np.sqrt( (state.zmax_r[j] - state.r_c[j])**2 + (state.zmax[j] - state.z_c[j])**2) #distance from the node to its corresponding normal

        state.nx_node[j] = state.nx_bt_bdy[j] * (1 - dist_bwd / (dist_bwd + dist_fwd) ) + state.nx_bt_bdy[j-1] * (dist_bwd / (dist_bwd + dist_fwd))
        state.nz_node[j] = state.nz_bt_bdy[j] * (1 - dist_bwd / (dist_bwd + dist_fwd) ) + state.nz_bt_bdy[j-1] * (dist_bwd / (dist_bwd + dist_fwd))

    norm = np.sqrt(state.nx_node**2 + state.nz_node**2)
    #normalize :
    state.nx_node = state.nx_node / norm
    state.nz_node = state.nz_node / norm

    state.tx_node = -state.nz_node
    state.tz_node = state.nx_node

def interpolate_normals(i, state, zM, bthy_m) :

    #bthy_m : bathy section of each ray (nray,)
    #zM : number of bounces nbounce < nray (nbounce,)
    #bthy_m[zM] : section where the bounce occurs (nbounce,)
    #nx_bt_bdy[bthy[zM]] : normal of the section where the bounce occurs (nbounce,)

    #distance between the two nodes
    dist_nodes = np.sqrt((state.zmax_r[bthy_m[zM]] - state.zmax_r[bthy_m[zM]+1])**2 + (state.zmax[bthy_m[zM]] - state.zmax[bthy_m[zM]+1])**2)

    #distance between the hitting point and the nodes
    dist_m = np.sqrt((state.r[i+1,zM] - state.zmax_r[bthy_m[zM]])**2+(state.z[i+1,zM] - state.zmax[bthy_m[zM]])**2)  #(zM,)
    dist_p =  np.sqrt((state.r[i+1,zM] - state.zmax_r[bthy_m[zM]+1])**2+(state.z[i+1,zM] - state.zmax[bthy_m[zM]+1])**2)  #(zM,)

    nx_bt_bdy = state.nx_node[bthy_m[zM]] * ( dist_p / dist_nodes) + state.nx_node[bthy_m[zM]+1] * (1- dist_p  / dist_nodes)
    nz_bt_bdy = state.nz_node[bthy_m[zM]] * ( dist_p / dist_nodes) + state.nz_node[bthy_m[zM]+1] * (1- dist_p  / dist_nodes)

    norm = np.sqrt(nx_bt_bdy**2 + nz_bt_bdy**2)

    nx_bt_bdy = nx_bt_bdy/norm
    nz_bt_bdy = nz_bt_bdy/norm

    state.ray_x_bdy[i+1,zM] = nx_bt_bdy
    state.ray_z_bdy[i+1,zM] = nz_bt_bdy

    return nx_bt_bdy, nz_bt_bdy


def crossing_depth(x1, y1, x2, y2, x3, y3, x4, y4) :
    z_bot = ((x1*y2 - y1*x2) * (y3 - y4) - (y1 - y2) * (x3*y4 - y3*x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    return z_bot
