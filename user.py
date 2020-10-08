"""File used to interact with the program"""

import core
import profile

#available exp : "TL", "R", "A"

def run() :
    for i in range(1) :


        print("Run " + str(i))
        print("####################")

        nr = 40

        param = {'nr' : nr, 'z0' : 0, 'zmin' : 0, 'rmin' : 0, 'rmax' : 20e3, 'zmax' : 5000,
                 'r0' : 0, 'angles' : (-60,-10), 'ds0' : 100, 'f' : 100, 'Lr' : 200, 'Lz' : 200,
                 'exp' : "R",
                 'r_rcvr' : 40e3, 'z_rcvr' : 250,
                 'compare_Bellhop' : 0, #Needs access to Bellhop, set it to False (0) if you trust this code
                 'speed_rand' : 0,
                 'speed_dist' : 'gaussian',
                 'speed_mean' : 0,
                 'speed_std' : 4,
                 'speed_Lz' : 100,
                 'speed_dim' : 1,
                 #'savefile' : "mcrun_"+str(i),
                 'savefile' : 'FLUOR_TL_0',
                 'plot' : 1, #When looping, should be set to 0
                 'save' : 0,
                 'exp_ind' : i,
                 'mean_prof' : 0,
                 'load_c' : 0,
                 'use_fortran' : 0, #speedup some functions, not ready yet
                 'range_dependent_bathy' : 1, #needs to load the bathy somewhere
                 'bathy_linterp' : 0 #interpolates linearly the bathymetry normals
        }

        #Rays unperturbed, ds0 = 200 : OK


        rays = core.core(param)


        rays.run()

if __name__ == '__main__' :
    run()
