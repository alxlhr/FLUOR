import numpy as np


def fluid_fluid(state,i,rind, theta_I) :

    #Snell's law of refraction :
        #k1 * cos(theta1) = k2 * cos(theta2)
        #theta1 = 90 - incident
        #k1 = omega/c1
        #k2 = omega/c2

        c1 = state.C[i+1,rind]
        c2 = state.c_bot

        k1 = 2*np.pi * state.f / c1
        k2 = 2*np.pi * state.f / c2

        theta1 = np.pi/2 - theta_I
        theta2 = np.arccos(k1/k2 * np.cos(theta1))

        #first approx :
        rho_1 = 1028 #kg/m**3


        #thetac = np.arccos(c1/c2)
        #theta_int = np.arctan(np.sqrt( (1 - (c2/c1)**2)/((state.rho_bot*c2/(rho_1*c1))**2 - 1)  ))#intromission angle

        Z1 = rho_1*c1 / np.sin(theta1)
        Z2 = state.rho_bot*c2 / np.sin(theta2)

        R = (Z2 - Z1) / (Z2 + Z1)
        T = 2*Z2 / (Z2 + Z1)

        #R[theta1 < theta_int] = 0
        #T[theta1 < theta_int] = 1
        #R = 1
        return R
