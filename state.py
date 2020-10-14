import numpy as np


class state(object) :
    """
    class containing the state of the variables
    """

    def __init__(self, param) :

        self.z0 = param['z0']
        self.r0 = param['r0']

        self.f = param['f']

        self.ds0 = param['ds0']
        self.nr = param['nr']

        #study area parameters
        self.zmax = param['zmax']
        self.zmin = param['zmin']
        self.rmax = param['rmax']
        self.rmin = param['rmin']
        self.Lr = param['Lr']
        self.Lz = param['Lz']

        self.angle_0 = np.deg2rad(np.linspace(param['angles'][0],param['angles'][1],self.nr))

        state.exp = param['exp']

        if state.exp == 'A' :
            self.r_rcvr = param['r_rcvr']
            self.z_rcvr = param['z_rcvr']

        self.n_max = np.int(np.ceil(param['rmax']/param['ds0'] * 1.5))

        #speed params
        self.speed_rand = param['speed_rand']
        self.speed_dist = param['speed_dist']
        if self.speed_dist != 'gaussian' :
            raise NameError('Wrong speed distribution')
        self.speed_mean = param['speed_mean']
        self.speed_std = param['speed_std']
        self.speed_Lz = param['speed_Lz']

        self.savefile = param['savefile']

        self.plot = param['plot']
        self.save = param['save']

        self.mean_prof = param['mean_prof']
        self.exp_ind = param['exp_ind']

        self.load_c = param['load_c']

        self.s_dim = param['speed_dim']

        self.use_fortran = param['use_fortran']

        self.rd_bathy = param['range_dependent_bathy']
        self.bathy_linterp = param['bathy_linterp']
        self.bathy_file = param['bathy_file']

        self.c_bot = param['c_bot']
        self.rho_bot = param['rho_bot']
