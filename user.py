"""File used to interact with the program"""

import core

#available exp : "TL", "R", "A"

for i in range(1) :

#i = 0

    print("Run " + str(i))
    print("####################")

    #nr = 25000
    nr = 1000

    param = {'nr' : nr, 'z0' : 1000, 'zmin' : 0, 'rmin' : 0, 'rmax' : 50000, 'zmax' : 5000,
    'r0' : 0, 'angles' : (-20,20), 'ds0' : 1000, 'f' : 1000, 'Lr' : 200, 'Lz' : 200,
    'exp' : "A",
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
    'load_c' : 0
    }

    #Rays unperturbed, ds0 = 200 : OK


    rays = core.core(param)


    rays.run()
