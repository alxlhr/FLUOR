import numpy as np
import loop
import speed
import bdy_loss
import matplotlib.pyplot as plt

def apply(i,state) :

    if state.rd_bathy == 0 :
        z_bot = state.zmax
        zM = (state.z[i+1,:] > z_bot)

        nz_bt_bdy = -1 * np.ones_like(state.z[i+1,:])
        tx_bt_bdy = -1 * np.ones_like(state.z[i+1,:])
        nx_bt_bdy = np.zeros_like(state.z[i+1,:])
        tz_bt_bdy = np.zeros_like(state.z[i+1,:])

    elif state.rd_bathy == 1 :

        #bthy_m = np.argmin(state.zmax_r[:,None] < state.r[i+1,:],axis = 0)-1

        bt_mp = np.argmin(state.zmax_r[:,None] < state.r[i+1,:],axis = 0)-1
        bt_mm =  np.argmin(state.zmax_r[:,None] < state.r[i,:],axis = 0)-1
        bt_mm[bt_mm < 0] = 0


        if np.all(bt_mm == bt_mp) :
            bthy_m = bt_mp.copy()

        #calculate interp point
        #check both segments, if it intercepts none, it doesn't cross the boundary
        else :
            [rb_m,zb_m] = crossing_depth(state.r[i,:],state.z[i,:],state.r[i+1,:],state.z[i+1,:],
            state.zmax_r[bt_mm], state.zmax[bt_mm], state.zmax_r[bt_mp], state.zmax[bt_mp])
            [rb_p,zb_p] = crossing_depth(state.r[i,:],state.z[i,:],state.r[i+1,:],state.z[i+1,:],
            state.zmax_r[bt_mp], state.zmax[bt_mp], state.zmax_r[bt_mp+1], state.zmax[bt_mp+1])

            dist_r_m = rb_m - state.zmax_r[bt_mm]
            dist_z_m = zb_m - state.zmax[bt_mm]
            dray_m = dist_r_m*state.tx_bt_bdy[bt_mm] + dist_z_m*state.tz_bt_bdy[bt_mm] #distance along the ray
            alpha_m = dray_m/state.bdy_dl[bt_mm]

            dist_r_p = rb_p - state.zmax_r[bt_mp]
            dist_z_p = zb_p - state.zmax[bt_mp]
            dray_p = dist_r_p*state.tx_bt_bdy[bt_mp] + dist_z_p*state.tz_bt_bdy[bt_mp] #distance along the ray
            alpha_p = dray_p/state.bdy_dl[bt_mp]

            bthy_m = bt_mp.copy()
            bthy_m[(alpha_m > 0) & (alpha_m < 1)] = bt_mm[(alpha_m > 0) & (alpha_m < 1)]

        bthy_M = bthy_m + 1

        #depth at the point (interpolate between bthy_m and bthy_M)
        B = np.abs(state.zmax_r[bthy_M] - state.zmax_r[bthy_m])

        #linear interpolation of the depth at that point (linterp from wikipedia)
        state.z_bot[i+1,:] = (state.zmax[bthy_m] * ( state.zmax_r[bthy_M] - state.r[i+1,:] ) + state.zmax[bthy_M] * ( state.r[i+1,:] - state.zmax_r[bthy_m] ) )/ B
        state.r_bot[i+1:] = state.zmax_r[bthy_M] - (state.zmax_r[bthy_M] - state.zmax_r[bthy_m]) * (state.z_bot[i+1,:] - state.zmax[bthy_M]) / (state.zmax[bthy_m] - state.zmax[bthy_M])
        z_bot = state.z_bot[i+1,:].copy()

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
            ds = state.ds0.copy()
            ds[zM] = (z_bot - state.z[i,zM]) / (state.C[i,zM]*state.Y[i,zM])
            loop.ray_step(i,zM,ds,state)
            nx_bt_bdy = nx_bt_bdy[zM]
            nz_bt_bdy = nz_bt_bdy[zM]
            tx_bt_bdy = tx_bt_bdy[zM]
            tz_bt_bdy = tz_bt_bdy[zM]
            kappa = np.zeros_like(nz_bt_bdy)
            #print(ds)
        elif state.rd_bathy == 1 :
            if state.bathy_linterp == 1 :
                reflexion_interp(state,zM,i,bthy_m, bthy_M)
                [nx_bt_bdy, nz_bt_bdy,tx_bt_bdy, tz_bt_bdy, kappa] = interpolate_normals(i, state,zM, bthy_m)
            elif state.bathy_linterp == 0 :
                reflexion_interp(state,zM,i,bthy_m, bthy_M)
                nx_bt_bdy = state.nx_bt_bdy[bthy_m[zM]]
                nz_bt_bdy = state.nz_bt_bdy[bthy_m[zM]]
                tx_bt_bdy = state.tx_bt_bdy[bthy_m[zM]]
                tz_bt_bdy = state.tz_bt_bdy[bthy_m[zM]]
                kappa = np.zeros_like(nx_bt_bdy)
            else :
                raise ValueError('Bad bathy interpolation option')

        dcdr = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],0,1, state.s_dim)
        dcdz = speed.get_der(state.f_interp,state.z[i+1,zM],state.r[i+1,zM],1,0, state.s_dim)

        #tangent : c * (X,Y)
        #normal :  c * (-Y,X)

        #Incident
        alpha = state.tx[i+1,zM]*nx_bt_bdy + nz_bt_bdy*state.tz[i+1,zM] #t_ray * boundary normal
        beta = tx_bt_bdy*state.tx[i+1,zM] + tz_bt_bdy*state.tz[i+1,zM] #t_ray * boundary tangent

#        print(tz_bt_bdy)

        nz_prev = state.nz[i+1,zM]
        nx_prev = state.nx[i+1,zM]
        tz_prev = state.tz[i+1,zM]
        tx_prev = state.tx[i+1,zM]

        #print(cn)

        #Reflected
        state.Y[i+1,zM] = state.Y[i+1,zM] - 2 * alpha * nz_bt_bdy / state.C[i+1,zM]
        state.X[i+1,zM] = state.X[i+1,zM] - 2 * alpha * nx_bt_bdy / state.C[i+1,zM]
        #state.X[i+1,zM] = (- alpha * nx_bt_bdy + beta * tx_bt_bdy) / state.C[i+1,zM]
        #state.Y[i+1,zM] = (- alpha * nz_bt_bdy + beta * tz_bt_bdy) / state.C[i+1,zM]

        state.tx[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]
        state.tz[i+1,zM] = state.C[i+1,zM]*state.Y[i+1,zM]
        state.nx[i+1,zM] = -state.C[i+1,zM]*state.Y[i+1,zM]
        state.nz[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]

        cn = dcdz * (nz_prev - state.nz[i+1,zM]) + dcdr * (nx_prev - state.nx[i+1,zM])
        cs = dcdz * (tz_prev - state.tz[i+1,zM]) + dcdr * (tx_prev - state.tx[i+1,zM])

        #cn = dcdz * ( - state.nz[i+1,zM])

        M = beta/alpha
        N = 2*kappa / (alpha * state.C[i+1,zM]**2) + M * (4*cn - 2*M*cs)/state.C[i+1,zM]**2

        #if zM[6] == True :
        #    print(tz_bt_bdy)

        #Ray normal and tangent
        theta_I = state.tx[i+1,zM]*nx_bt_bdy + state.tz[i+1,zM]*nz_bt_bdy

        q_prev = state.q[i+1,zM].copy()
        state.q[i+1,zM] = q_prev
        state.p[i+1,zM] = state.p[i+1,zM] + q_prev * np.abs(N)


        #boundary losses
        R = bdy_loss.fluid_fluid(state,i,zM,theta_I)
        state.amp[i+1,zM] = state.amp[i,zM] * np.abs(R)
        state.phi[i+1,zM] = state.phi[i,zM] + np.angle(R)

        #Angles check

        theta_R = state.tx[i+1,zM]*nx_bt_bdy + state.tz[i+1,zM]*nz_bt_bdy
        #print("**********")
        #print("I : ", theta_I)
        #print("R : ", theta_R)
        chck_ang = theta_I - theta_R < 1e-12
        if (np.any(chck_ang == False)) :
            print("Problem with incident/reflected angles : ",i)

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

"""
def recalculate_step(state, i, zM, bthy_m, bthy_M, nx_bt_bdy, nz_bt_bdy) :

    d0_r = state.r[i,zM] - state.zmax_r[bthy_m[zM]]
    d0_z = state.z[i,zM] - state.zmax[bthy_m[zM]]

    d_r = state.r[i+1,zM] - state.zmax_r[bthy_m[zM]]
    d_z = state.z[i+1,zM] - state.zmax[bthy_m[zM]]

    de_0 = d0_r * nx_bt_bdy + d0_z * nz_bt_bdy
    de = d_r * nx_bt_bdy + d_z * nz_bt_bdy

    dh = de_0 / (de_0 + np.abs(de)) * state.ds0[zM]

    ds = np.zeros_like(state.ds0)
    ds[zM] = dh#de_0 / np.abs(state.tx[i+1,zM]*nx_bt_bdy + state.tz[i+1,zM]*nz_bt_bdy)

    return ds
"""
def reflexion_interp(state,zM,i,bthy_m, bthy_M) :

    [rb,zb] = crossing_depth(state.r[i,zM],state.z[i,zM],state.r[i+1,zM],state.z[i+1,zM], state.zmax_r[bthy_m[zM]], state.zmax[bthy_m[zM]], state.zmax_r[bthy_M[zM]], state.zmax[bthy_M[zM]])

    dist_r = rb - state.r[i,zM]
    dist_z = zb - state.z[i,zM]

    dray = np.sqrt(dist_r**2 + dist_z**2) #distance along the ray

    alpha = 1 - dray/state.ds0[zM] #linterp coef
    #alpha[alpha > 1] = 1
    #alpha = 1 - alpha

    state.r[i+1,zM] = rb
    state.z[i+1,zM] = zb

    state.X[i+1,zM] = state.X[i,zM]*alpha + (1-alpha)*state.X[i+1,zM]
    state.Y[i+1,zM] = state.Y[i,zM]*alpha + (1-alpha)*state.Y[i+1,zM]

    state.p[i+1,zM] = state.p[i,zM]*alpha + (1-alpha)*state.p[i+1,zM]
    state.q[i+1,zM] = state.q[i,zM]*alpha + (1-alpha)*state.q[i+1,zM]

    state.T[i+1,zM] = state.T[i,zM]*alpha + (1-alpha)*state.T[i+1,zM]
    state.angle[i+1,zM] = np.arctan((state.z[i+1,zM] - state.z[i,zM])/(state.r[i+1,zM] - state.r[i,zM]))
    state.angle[i+1,np.isnan(state.angle[i+1,:])] = np.pi/2
    state.C[i+1,zM] = speed.get_speed(state.z[i+1,zM], state.r[i+1,zM], state.f_interp, state.s_dim)

    state.tx[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]
    state.tz[i+1,zM] = state.C[i+1,zM]*state.Y[i+1,zM]
    state.nx[i+1,zM] = -state.C[i+1,zM]*state.Y[i+1,zM]
    state.nz[i+1,zM] = state.C[i+1,zM]*state.X[i+1,zM]

def calculate_normals(state) :

    dzmax = np.diff(state.zmax)
    drmax = np.diff(state.zmax_r)

    alpha_r = np.arctan(dzmax/drmax)

    #normal to the bathy section
    state.nx_bt_bdy = np.sin(alpha_r)
    state.nz_bt_bdy = -np.cos(alpha_r) #z > 0 upward

    state.tx_bt_bdy = state.nz_bt_bdy
    state.tz_bt_bdy = -state.nx_bt_bdy #(t,n,y) must form a right handed coordinate system, y going towards the screen

    #get center of the bathy sections

    state.bdy_dl = np.sqrt(dzmax**2 + drmax**2)
    state.z_c = state.bdy_dl/2 * np.sin(alpha_r) + state.zmax[:-1]
    state.r_c =  state.bdy_dl/2 * np.cos(alpha_r) + state.zmax_r[:-1]

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

    state.tx_node = state.nz_node
    state.tz_node = -state.nx_node

    normal_angle = np.arctan(state.nz_node/np.abs(state.nx_node))
    dist_nodes = np.sqrt((state.zmax_r[1:] - state.zmax_r[:-1])**2 + (state.zmax[1:] - state.zmax[:-1])**2)

    state.kappa = np.zeros_like(state.zmax_r)
    state.kappa[:-1] = np.diff(normal_angle)/dist_nodes

def interpolate_normals(i, state, zM, bthy_m) :

    #bthy_m : bathy section of each ray (nray,)
    #zM : number of bounces nbounce < nray (nbounce,)
    #bthy_m[zM] : section where the bounce occurs (nbounce,)
    #nx_bt_bdy[bthy[zM]] : normal of the section where the bounce occurs (nbounce,)
    #print(zM)
    #distance between the two nodes
    dist_nodes = np.sqrt((state.zmax_r[bthy_m[zM]] - state.zmax_r[bthy_m[zM]+1])**2 + (state.zmax[bthy_m[zM]] - state.zmax[bthy_m[zM]+1])**2)

    #distance between the hitting point and the nodes
    dist_m = np.sqrt((state.r[i+1,zM] - state.zmax_r[bthy_m[zM]])**2+(state.z[i+1,zM] - state.zmax[bthy_m[zM]])**2)  #(zM,)
    dist_p =  np.sqrt((state.r[i+1,zM] - state.zmax_r[bthy_m[zM]+1])**2+(state.z[i+1,zM] - state.zmax[bthy_m[zM]+1])**2)  #(zM,)

    #linear interpolation
    nx_bt_bdy = state.nx_node[bthy_m[zM]+1] * (1-dist_p / dist_nodes) + state.nx_node[bthy_m[zM]] * (dist_p  / dist_nodes)
    nz_bt_bdy = state.nz_node[bthy_m[zM]+1] * (1-dist_p / dist_nodes) + state.nz_node[bthy_m[zM]] * (dist_p  / dist_nodes)
    tx_bt_bdy = state.tx_node[bthy_m[zM]+1] * (1-dist_p / dist_nodes) + state.tx_node[bthy_m[zM]] * (dist_p  / dist_nodes)
    tz_bt_bdy = state.tz_node[bthy_m[zM]+1] * (1-dist_p / dist_nodes) + state.tz_node[bthy_m[zM]] * (dist_p  / dist_nodes)
    kappa = state.kappa[bthy_m[zM]+1] * (1-dist_p / dist_nodes) + state.kappa[bthy_m[zM]] * (dist_p  / dist_nodes)

    norm = np.sqrt(nx_bt_bdy**2 + nz_bt_bdy**2)

    nx_bt_bdy = nx_bt_bdy/norm
    nz_bt_bdy = nz_bt_bdy/norm

    state.ray_x_bdy[i+1,zM] = nx_bt_bdy
    state.ray_z_bdy[i+1,zM] = nz_bt_bdy

    return nx_bt_bdy, nz_bt_bdy, tx_bt_bdy, tz_bt_bdy, kappa


def crossing_depth(x1, y1, x2, y2, x3, y3, x4, y4) :
    "From wikipedia"
    z_bot = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    r_bot = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4) ) / ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4) )
    return r_bot,z_bot
