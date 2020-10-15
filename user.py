"""File used to interact with the program"""

import core
import profile
"""
from os import *
import sys
import numpy as np
from scipy.io import *
import matplotlib.pyplot as plt
sys.path.append ("/home/alexandre/Documents/Stage_M2/Orlando_Python/at/Python/")
from readshd import *
from plotray import *
"""
#available exp : "TL", "R", "A"

def run() :
    for i in range(1) :
        i = 3

        print("Run " + str(i))
        print("####################")

        nr = 20

        param = {'nr' : nr, 'z0' : 10, 'zmin' : 0, 'rmin' : 0, 'rmax' : 20e3, 'zmax' : 1.2e3,
                 'r0' : 0, 'angles' : (-20,20), 'ds0' : 20, 'f' : 1000, 'Lr' : 200, 'Lz' : 200,
                 'exp' : "R",
                 'r_rcvr' : 19e3, 'z_rcvr' : 100,
                 'compare_Bellhop' : 0, #Needs access to Bellhop, set it to False (0) if you trust this code
                 'speed_rand' : 0,
                 'speed_dist' : 'gaussian',
                 'speed_mean' : 0,
                 'speed_std' : 4,
                 'speed_Lz' : 100,
                 'speed_dim' : 1,
                 #'savefile' : "mcrun_"+str(i),
                 'savefile' : 'FLUOR_Bathy_'+str(i),
                 'plot' : 1, #When looping, should be set to 0
                 'save' : 0,
                 'exp_ind' : i,
                 'mean_prof' : 0,
                 'load_c' : 1,
                 'use_fortran' : 0, #speedup some functions, not ready yet
                 'range_dependent_bathy' : 1, #needs to load the bathy somewhere
                 'bathy_linterp' : 0, #interpolates linearly the bathymetry normals
                 'bathy_file' : '../gebco_2020_n38.3_s37.8_w15.3_e15.8.nc',
                 'c_bot' : 1600,
                 'rho_bot' : 2000
                 #'bathy_file' : '/home6/lheral/Documents/Stage_M2/gebco_2020_n38.3_s37.8_w15.3_e15.8.nc'
        }

        #Rays unperturbed, ds0 = 200 : OK


        rays = core.core(param)


        rays.run()


if __name__ == '__main__' :
    run()
