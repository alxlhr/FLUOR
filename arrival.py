import numpy as np

def calc_arr(state) :
    state.W = np.abs(state.q * state.dangle0[np.newaxis,:])

    state.Angle_rcvr = np.zeros((state.nr))
    state.Delay_rcvr = np.zeros((state.nr))
    state.Amp_rcvr = np.zeros((state.nr))
    #state.Amp_rays = np.zeros((state.nr), dtype = complex)
    state.ray_num = np.zeros((state.nr), dtype = np.int32)


    state.eigen_ray = 0

    loop(state)

    n_arr = state.Angle_rcvr[:state.eigen_ray]
    state.Angle_rcvr = n_arr

    n_arr = state.Delay_rcvr[:state.eigen_ray]
    state.Delay_rcvr = n_arr

    n_arr = state.Amp_rcvr[:state.eigen_ray]
    state.Amp_rcvr = n_arr

    n_arr = state.ray_num[:state.eigen_ray]
    state.ray_num = n_arr

    print('     Arrivals done       ')
    print(' ')


def loop(state) :

    for j in range(state.nr) :
        for i in range(1,state.n_max) :
            rec = ((state.r_rcvr >= state.r[i-1,j]) & (state.r_rcvr < state.r[i,j]))
            if rec == True :

                n_ = np.abs( (state.r_rcvr - state.r[i,j])*state.nx[i,j] + (state.z_rcvr - state.z[i,j])*state.nz[i,j])
                s_ = (state.r_rcvr - state.r[i,j])*state.tx[i,j] + (state.z_rcvr - state.z[i,j])*state.tz[i,j]
                n_val = (n_ <= 2.0)

                if n_val == True :

                    #print(state.q[i,j])

                    length_ray = np.sqrt((state.r[i,j] - state.r[i-1,j])**2 + (state.z[i,j] - state.z[i-1,j])**2)
                    al = s_ / length_ray

                    T_ = state.T[i-1,j] + al * (state.T[i,j] - state.T[i-1,j])
                    q_ = state.q[i-1,j] + al * (state.q[i,j] - state.q[i-1,j])
                    W_ = state.W[i-1,j] + al * (state.W[i,j] - state.W[i-1,j])
                    r_ = state.r[i-1,j] + al * (state.r[i,j] - state.r[i-1,j])
                    angle_ = state.angle[i-1,j] + al * (state.angle[i,j] - state.angle[i-1,j])

                    A_ = (-1j)**state.m[i,j] * np.sqrt(np.abs(state.C[i,j]*np.cos(state.angle_0[j])/(r_*state.c0*q_)))

                    state.Angle_rcvr[state.eigen_ray] = np.rad2deg(angle_)
                    state.Delay_rcvr[state.eigen_ray] = T_

                    state.Amp_rcvr[state.eigen_ray] = np.abs(A_)#np.abs((W_ - n_)/W_*A_)#*np.exp(1j*2*np.pi*f*T_))
                    #state.Amp_rays[state.eigen_ray] = (W_ - n_)/W_*A_*np.exp(1j*2*np.pi*state.f*T_)

                    state.ray_num[state.eigen_ray] = j

                    state.eigen_ray = state.eigen_ray + 1
