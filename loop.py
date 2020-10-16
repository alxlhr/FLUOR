import numpy as np
import speed


def ray_step(i,j,ds,state) :

    #print(np.shape(ds))
    out = (state.r[i,:] < state.rmax)
    #out = (state.r[i,:] > state.rmin) & out 
    #print(out)
    j = j & out
    ds = ds[j]

    """f(tn,yn)"""
    state.C[i,j] = speed.get_speed(state.z[i,j], state.r[i,j], state.f_interp, state.s_dim)

    dr_1 = state.C[i,j]*state.X[i,j]
    dz_1 = state.C[i,j]*state.Y[i,j]

    dcdr_1 = speed.get_der(state.f_interp,state.z[i,j],state.r[i,j],0,1, state.s_dim)
    dcdz_1 = speed.get_der(state.f_interp,state.z[i,j],state.r[i,j],1,0, state.s_dim)

    #k1
    dX_1 = - 1/state.C[i,j]**2 * dcdr_1
    dY_1 = - 1/state.C[i,j]**2 * dcdz_1

    """f(yn+h/2*k1)"""
    #y_n + h/2 * k1
    r_2 = state.r[i,j] + ds/2 * dr_1
    z_2 = state.z[i,j] + ds/2 * dz_1

    C_2 = speed.get_speed(z_2,r_2, state.f_interp, state.s_dim)

    Y_2 = state.Y[i,j] + ds/2 * dY_1
    X_2 = state.X[i,j] + ds/2 * dX_1

    dr_2 = C_2*X_2
    dz_2 = C_2*Y_2

    dcdr_2 = speed.get_der(state.f_interp,z_2,r_2,0,1, state.s_dim)
    dcdz_2 = speed.get_der(state.f_interp,z_2,r_2,1,0, state.s_dim)

    #k2
    dX_2 = - 1/C_2**2 * dcdr_2
    dY_2 = - 1/C_2**2 * dcdz_2

    """f(yn+h/2*k2)"""
    #y_n + h/2 * k1
    r_3 = state.r[i,j] + ds/2 * dr_2
    z_3 = state.z[i,j] + ds/2 * dz_2

    C_3 = speed.get_speed(z_3,r_3, state.f_interp, state.s_dim)

    Y_3 = state.Y[i,j] + ds/2 * dY_2
    X_3 = state.X[i,j] + ds/2 * dX_2

    dr_3 = C_3*X_3
    dz_3 = C_3*Y_3

    dcdr_3 = speed.get_der(state.f_interp,z_3,r_3,0,1, state.s_dim)
    dcdz_3 = speed.get_der(state.f_interp,z_3,r_3,1,0, state.s_dim)

    #print(dcdz_3)

    #k3
    dX_3 = - 1/C_3**2 * dcdr_3
    dY_3 = - 1/C_3**2 * dcdz_3

    """f(yn+h*k3)"""
    #y_n + h/2 * k1
    r_4 = state.r[i,j] + ds * dr_3
    z_4 = state.z[i,j] + ds * dz_3

    C_4 = speed.get_speed(z_4,r_4, state.f_interp, state.s_dim)

    Y_4 = state.Y[i,j] + ds * dY_3
    X_4 = state.X[i,j] + ds * dX_3

    dr_4 = C_4*X_4
    dz_4 = C_4*Y_4

    dcdr_4 = speed.get_der(state.f_interp,z_4,r_4,0,1, state.s_dim)
    dcdz_4 = speed.get_der(state.f_interp,z_4,r_4,1,0, state.s_dim)

    #k4
    dX_4 = - 1/C_4**2 * dcdr_4
    dY_4 = - 1/C_4**2 * dcdz_4

    #i+1
    state.r[i+1,j] = state.r[i,j] + ds/6 * (dr_1 + 2*dr_2 + 2*dr_3 + dr_4)
    state.z[i+1,j] = state.z[i,j] + ds/6 * (dz_1 + 2*dz_2 + 2*dz_3 + dz_4)

    state.X[i+1,j] = state.X[i,j] + ds/6* (dX_1 + 2*dX_2 + 2*dX_3 + dX_4)
    state.Y[i+1,j] = state.Y[i,j] + ds/6* (dY_1 + 2*dY_2 + 2*dY_3 + dY_4)  #zeta

    #TODO : bad practice, switch to masked array (worth it ?)
    state.r[i+1,np.invert(out)] = np.nan
    state.z[i+1,np.invert(out)] = np.nan

    state.X[i+1,np.invert(out)] = np.nan
    state.Y[i+1,np.invert(out)] = np.nan  #zeta

    #Dynamic rays

    """f(tn,yn)"""
    dcdr2_1 = speed.get_der(state.f_interp,state.z[i,j],state.r[i,j],0,2, state.s_dim)
    dcdz2_1 = speed.get_der(state.f_interp,state.z[i,j],state.r[i,j],2,0, state.s_dim)
    dcdrdz_1 = speed.get_der(state.f_interp,state.z[i,j],state.r[i,j],1,1, state.s_dim)

    cnn_1 = (dcdz2_1* state.X[i,j]**2 + dcdr2_1 * state.Y[i,j]**2 - 2*dcdrdz_1*state.X[i,j]*state.Y[i,j])

    dp_1 =  - cnn_1 * state.q[i,j]
    dq_1 = state.C[i,j]*state.p[i,j]

    """k2"""

    p_2 = state.p[i,j] + dp_1*ds/2
    q_2 = state.q[i,j] + dq_1*ds/2

    dcdr2_2 = speed.get_der(state.f_interp,z_2,r_2,0,2, state.s_dim)
    dcdz2_2 = speed.get_der(state.f_interp,z_2,r_2,2,0, state.s_dim)
    dcdrdz_2 = speed.get_der(state.f_interp,z_2,r_2,1,1, state.s_dim)

    cnn_2 = (dcdz2_2* X_2**2 + dcdr2_2 * Y_2**2 - 2*dcdrdz_2*X_2*Y_2)

    dp_2 =  - cnn_2 * q_2
    dq_2 = C_2*p_2

    """k3"""

    p_3 = p_2 + dp_2*ds/2
    q_3 = q_2 + dq_2*ds/2

    dcdr2_3 = speed.get_der(state.f_interp,z_3,r_3,0,2, state.s_dim)
    dcdz2_3 = speed.get_der(state.f_interp,z_3,r_3,2,0, state.s_dim)
    dcdrdz_3 = speed.get_der(state.f_interp,z_3,r_3,1,1, state.s_dim)

    cnn_3 = (dcdz2_3* X_3**2 + dcdr2_3 * Y_3**2 - 2*dcdrdz_3*X_3*Y_3)

    dp_3 =  - cnn_3 * q_3
    dq_3 = C_3*p_3

    """k4"""

    p_4 = p_3 + dp_3*ds
    q_4 = q_3 + dq_3*ds

    dcdr2_4 = speed.get_der(state.f_interp,z_4,r_4,0,2, state.s_dim)
    dcdz2_4 = speed.get_der(state.f_interp,z_4,r_4,2,0, state.s_dim)
    dcdrdz_4 = speed.get_der(state.f_interp,z_4,r_4,1,1, state.s_dim)

    cnn_4 = (dcdz2_4* X_4**2 + dcdr2_4 * Y_4**2 - 2*dcdrdz_4*X_4*Y_4)

    dp_4 =  - cnn_4 * q_4
    dq_4 = C_4*p_4


    state.p[i+1,j] = state.p[i,j] + ds/6 * (dp_1 + 2*dp_2 + 2*dp_3 + dp_4)
    state.q[i+1,j] = state.q[i,j] + ds/6 * (dq_1 + 2*dq_2 + 2*dq_3 + dq_4)

    #Phase

    dT_1 = 1/state.C[i,j]
    dT_2 = 1/C_2
    dT_3 = 1/C_3
    dT_4 = 1/C_4

    state.T[i+1,j] = state.T[i,j] + ds/6 * (dT_1 + 2*dT_2 + 2*dT_3 + dT_4)

    #print('i : ', i)
    #print('r : ',state.r[i+1,j])
    state.angle[i+1,j] = np.arctan((state.z[i+1,j] - state.z[i,j])/(state.r[i+1,j] - state.r[i,j]))
    #print('angle_ : ', (state.r[i+1,j] - state.r[i,j]))

    state.C[i+1,j] = speed.get_speed(state.z[i+1,j], state.r[i+1,j], state.f_interp, state.s_dim)

    state.tx[i+1,j] = state.C[i+1,j]*state.X[i+1,j]
    state.tz[i+1,j] = state.C[i+1,j]*state.Y[i+1,j]
    state.nx[i+1,j] = -state.C[i+1,j]*state.Y[i+1,j]
    state.nz[i+1,j] = state.C[i+1,j]*state.X[i+1,j]

    state.amp[i+1,j] = state.amp[i,j]
    state.phi[i+1,j] = state.phi[i,j]

    #print(state.tz[i+1,j])

    """
    print("*****")
    print(i)
    print(state.tx[i+1,j])
    
    """
