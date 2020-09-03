import netCDF4 as nc
import numpy as np
import os.path as path
import os

def save(state, params) :
    if state.speed_rand == 1 or state.speed_mean == 1:
        dir = '/media/alexandre/DATA/Stage_M2/exp/arrivals/dist_z/2D_arr/'
    else :
        #dir = '/media/alexandre/DATA/Stage_M2/exp/arrivals/dist_z/rec_c_500_100_n3/'
        dir = '/media/alexandre/DATA/Stage_M2/exp/arrivals/dist_z/2D_arr/'

    if (path.exists(dir) == False) :
        os.mkdir(dir)

    if (path.exists(dir + state.savefile+'.nc')) :
        print('# WARNING: savefile already exists, changing its name')
        state.savefile = state.savefile+'_bis'

    ncout = nc.Dataset(dir+state.savefile+'.nc','w')

    n_rays = ncout.createDimension('n_rays', state.nr)
    n_points = ncout.createDimension('n_points', state.n_max)

    Lr = ncout.createDimension('Lr', state.Lr)
    Lz = ncout.createDimension('Lz', state.Lz)

    #parameters
    parms = ncout.createGroup('params')
    for k,v in params.items() :
        setattr(parms, k, v)

    #Celerities
    if state.s_dim == 1 :
        C_bkgnd = ncout.createVariable('C_bkgnd', state.c_bk.dtype, dimensions=('Lz'))
    elif state.s_dim == 2 :
        C_bkgnd = ncout.createVariable('C_bkgnd', state.c_bk.dtype, dimensions=('Lz','Lr'))
    C_rays = ncout.createVariable('C_rays', state.C.dtype, dimensions=('n_points','n_rays'))
    print(np.shape(C_bkgnd))
    print(np.shape(state.c_bk))
    C_bkgnd[...] = state.c_bk
    C_rays[...] = state.C

    #Ray variables

    if state.exp != 'A' :

        X = ncout.createVariable('X', state.X.dtype, dimensions=('n_points','n_rays'))
        Y = ncout.createVariable('Y', state.Y.dtype, dimensions=('n_points','n_rays'))
        r = ncout.createVariable('r', state.r.dtype, dimensions=('n_points','n_rays'))
        z = ncout.createVariable('z', state.z.dtype, dimensions=('n_points','n_rays'))
        p = ncout.createVariable('p', state.p.dtype, dimensions=('n_points','n_rays'))
        q = ncout.createVariable('q', state.q.dtype, dimensions=('n_points','n_rays'))
        T = ncout.createVariable('T', state.T.dtype, dimensions=('n_points','n_rays'))
        W = ncout.createVariable('W', state.W.dtype, dimensions=('n_points','n_rays'))
        angle = ncout.createVariable('angle', state.angle.dtype, dimensions=('n_points','n_rays'))
        m = ncout.createVariable('m', state.m.dtype, dimensions=('n_points','n_rays'))
        tx = ncout.createVariable('tx', state.tx.dtype, dimensions=('n_points','n_rays'))
        tz = ncout.createVariable('tz', state.tz.dtype, dimensions=('n_points','n_rays'))
        nx = ncout.createVariable('nx', state.nx.dtype, dimensions=('n_points','n_rays'))
        nz = ncout.createVariable('nz', state.nz.dtype, dimensions=('n_points','n_rays'))
        bdy_bot = ncout.createVariable('bdy_bot', state.nz.dtype, dimensions=('n_rays'))
        bdy_top = ncout.createVariable('bdy_top', state.nz.dtype, dimensions=('n_rays'))

        X[...] = state.X
        Y[...] = state.Y
        r[...] = state.r
        z[...] = state.z
        p[...] = state.p
        q[...] = state.q
        T[...] = state.T
        W[...] = state.W
        angle[...] = state.angle
        m[...] = state.m
        tx[...] = state.tx
        tz[...] = state.tz
        nx[...] = state.nx
        nz[...] = state.nz
        bdy_bot[...] = state.bdy_bot
        bdy_top[...] = state.bdy_top

    #TL
    if state.exp == 'TL' :
        P_real = ncout.createVariable('P_real', state.P.real.dtype, dimensions=('Lz','Lr'))
        P_imag = ncout.createVariable('P_imag', state.P.imag.dtype, dimensions=('Lz','Lr'))
        TL = ncout.createVariable('TL', state.TL.dtype, dimensions=('Lz','Lr'))

        P_real[...] = state.P.real
        P_imag[...] = state.P.imag
        TL[...] = state.TL

    #A
    if state.exp == 'A' :

        n_eigen = ncout.createDimension('n_eigen', state.eigen_ray)

        Angle_rcvr = ncout.createVariable('Angle_rcvr', state.Angle_rcvr.dtype, dimensions=('n_eigen',))
        Amp_rcvr = ncout.createVariable('Amp_rcvr', state.Amp_rcvr.dtype, dimensions=('n_eigen',))
        #Amp_rays = ncout.createVariable('Amp_rays', state.Amp_rays.dtype, dimensions=('n_eigen'))
        ray_ids = ncout.createVariable('ray_ids', state.ray_num.dtype, dimensions=('n_eigen',))
        Delay_rcvr = ncout.createVariable('Delay_rcvr', state.Delay_rcvr.dtype, dimensions = ('n_eigen',))

        Angle_rcvr[...] = state.Angle_rcvr
        Amp_rcvr[...] = state.Amp_rcvr
        #Amp_rays[...] = state.Amp_rays
        ray_ids[...] = state.ray_num
        Delay_rcvr[...] = state.Delay_rcvr

        X = ncout.createVariable('X', state.X.dtype, dimensions=('n_points','n_eigen'))
        Y = ncout.createVariable('Y', state.Y.dtype, dimensions=('n_points','n_eigen'))
        r = ncout.createVariable('r', state.r.dtype, dimensions=('n_points','n_eigen'))
        z = ncout.createVariable('z', state.z.dtype, dimensions=('n_points','n_eigen'))
        p = ncout.createVariable('p', state.p.dtype, dimensions=('n_points','n_eigen'))
        q = ncout.createVariable('q', state.q.dtype, dimensions=('n_points','n_eigen'))
        T = ncout.createVariable('T', state.T.dtype, dimensions=('n_points','n_eigen'))
        W = ncout.createVariable('W', state.W.dtype, dimensions=('n_points','n_eigen'))
        angle = ncout.createVariable('angle', state.angle.dtype, dimensions=('n_points','n_eigen'))
        m = ncout.createVariable('m', state.m.dtype, dimensions=('n_points','n_eigen'))
        tx = ncout.createVariable('tx', state.tx.dtype, dimensions=('n_points','n_eigen'))
        tz = ncout.createVariable('tz', state.tz.dtype, dimensions=('n_points','n_eigen'))
        nx = ncout.createVariable('nx', state.nx.dtype, dimensions=('n_points','n_eigen'))
        nz = ncout.createVariable('nz', state.nz.dtype, dimensions=('n_points','n_eigen'))
        bdy_bot = ncout.createVariable('bdy_bot', state.bdy_bot.dtype, dimensions=('n_eigen'))
        bdy_top = ncout.createVariable('bdy_top', state.bdy_top.dtype, dimensions=('n_eigen'))

        X[...] = state.X[:,state.ray_num]
        Y[...] = state.Y[:,state.ray_num]
        r[...] = state.r[:,state.ray_num]
        z[...] = state.z[:,state.ray_num]
        p[...] = state.p[:,state.ray_num]
        q[...] = state.q[:,state.ray_num]
        T[...] = state.T[:,state.ray_num]
        W[...] = state.W[:,state.ray_num]
        angle[...] = state.angle[:,state.ray_num]
        m[...] = state.m[:,state.ray_num]
        tx[...] = state.tx[:,state.ray_num]
        tz[...] = state.tz[:,state.ray_num]
        nx[...] = state.nx[:,state.ray_num]
        nz[...] = state.nz[:,state.ray_num]
        bdy_bot[...] = state.bdy_bot[state.ray_num]
        bdy_top[...] = state.bdy_top[state.ray_num]
