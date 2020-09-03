import numpy as np


def step(i,state) :

    if (i > 1) :
        dsign = np.sign(state.q[i-1,:]) - np.sign(state.q[i,:])
        s_ind = (dsign != 0)
        if np.size(state.m[i,s_ind]) :
            state.m[i,:] = state.m[i-1,:]
            state.m[i,s_ind] = state.m[i-1,s_ind] + 1
        else :
            state.m[i,:] = state.m[i-1,:]
