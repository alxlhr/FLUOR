import numpy as np
import state
import speed
import loop
import caustics
import plot
import boundary
import loss
import timeit
import arrival
import IO
import gc
import netCDF4 as nc
import matplotlib.pyplot as plt

class core(object) :

    "core class of the ray code"

    def __init__(self, param) :

        self.state = state.state(param)
        self.params = param
        self.init_speed()
        self.init_var(param)
        self.check_res()


    def init_speed(self) :

        z_bk = np.linspace(self.state.zmin-10,self.state.zmax+10,self.state.Lz)
        r_bk = np.linspace(self.state.rmin,self.state.rmax,self.state.Lr)

        if self.state.load_c == 1 :
            ncf=nc.Dataset('/home/alexandre/Documents/Stage_M2/Python/EOF/rec_cel_m1.nc')
            z_int = ncf.variables['z_rec'][...]
            ncf.close()

        else :
            z_int = z_bk

        z_bk = np.linspace(self.state.zmin-10,self.state.zmax+10,self.state.Lz)
        r_bk = np.linspace(self.state.rmin,self.state.rmax,self.state.Lr)
        self.state.c_bk = np.squeeze(speed.GenerateRandom(z_int,r_bk,self.state))
        self.state.f_interp = speed.Init(z_int,r_bk,self.state.c_bk, z_bk, self.state.s_dim)

    def init_var(self,param) :
        self.state.c0 = speed.get_speed(param['z0'],param['r0'], self.state.f_interp, self.state.s_dim)

        X0 = np.cos(self.state.angle_0)/self.state.c0
        Y0 = np.sin(self.state.angle_0)/self.state.c0

        #Rays
        self.state.X = np.zeros((self.state.n_max,self.state.nr))
        self.state.Y = np.zeros((self.state.n_max,self.state.nr))
        self.state.r = np.zeros((self.state.n_max,self.state.nr))
        self.state.z = np.zeros((self.state.n_max,self.state.nr))
        self.state.C = np.zeros((self.state.n_max,self.state.nr))
        self.state.bdy_top = np.zeros((self.state.nr))
        self.state.bdy_bot = np.zeros((self.state.nr))

        #TL
        self.state.p = np.zeros((self.state.n_max,self.state.nr))
        self.state.q = np.zeros((self.state.n_max,self.state.nr))
        self.state.T = np.zeros((self.state.n_max,self.state.nr))
        self.state.W = np.zeros((self.state.n_max,self.state.nr))
        self.state.angle = np.zeros((self.state.n_max,self.state.nr))
        self.state.m = np.zeros((self.state.n_max,self.state.nr))
        self.state.dangle0 = np.zeros((self.state.nr))

        self.state.X[0,:] = X0
        self.state.Y[0,:] = Y0
        self.state.r[0,:] = self.state.r0
        self.state.z[0,:] = self.state.z0
        self.state.p[0,:] = 1/self.state.c0
        self.state.T[0,:] = 0
        self.state.C[0,:] = self.state.c0
        self.state.angle[0,:] = self.state.angle_0
        self.state.dangle0[:-1] = self.state.angle[0,1:] - self.state.angle[0,:-1]

    def check_res(self) :
        temp = self.state.dangle0 * self.state.rmax / np.cos(np.max(np.abs(self.state.angle_0)))
        hz = (self.state.zmax - self.state.zmin)/self.state.Lz

        if np.max(temp) >  hz :
            print('delta z is too big : increase vertical resolution or the number of rays')




    def run(self) :
        starttime = timeit.default_timer()

        print(self.state.exp)
        if self.state.exp == "R" :
            print("### Calculating Rays ###")
        elif self.state.exp == "TL" :
            print("### Calculating TL ###")
        elif self.state.exp == "A" :
            print("### Calculating Arrivals ###")
        else :
            raise NameError('Wrong exp name')

        def_arr = np.ones_like(self.state.z[0,:],dtype = bool)
        for i in range(self.state.n_max-1) :
            if (i%10 == 0) :
                print("%i / %i" %(i, self.state.n_max), end = '\r')
            loop.ray_step(i,def_arr,self.state.ds0,self.state)
            caustics.step(i,self.state)
            boundary.apply(i,self.state)
            gc.collect()

        loss.calc_normals(self.state)
        if self.state.exp == "TL" :
            loss.calc_TL(self.state)
        if self.state.exp == "A" :
            arrival.calc_arr(self.state)

        print("Elapsed time : %i sec" % (timeit.default_timer() - starttime))

        if self.state.save == 1 :
            IO.save(self.state, self.params)

        if self.state.plot == 1 :
            plot.show(self.state)

        gc.collect()
