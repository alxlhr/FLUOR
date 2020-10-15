import numpy as np
import scipy.interpolate as interp
import netCDF4 as nc
import matplotlib.pyplot as plt



def load(state) :
    #Load the bathymetry 
    nc_bathy = nc.Dataset(state.bathy_file)
    lat = nc_bathy.variables['lat'][...]
    lon = nc_bathy.variables['lon'][...]
    bathy = nc_bathy.variables['elevation'][...]

    Lon,Lat = np.meshgrid(lon,lat)
    L0 = np.argmin(np.abs(Lat - 38.05), axis = 0)
    state.zmax = -1*bathy[L0][0,:]
    indz = state.zmax > 0
    state.zmax = state.zmax[indz]
    state.zmax_r = lon[indz]*111e3
    state.zmax_r = state.zmax_r - np.min(state.zmax_r)
   
    state.rmin = np.min(state.zmax_r)
    state.rmax = np.max(state.zmax_r)
#state.zmax_r = np.linspace(state.rmin,state.rmax,len(state.zmax))


def interpolate(state,r_int) :
    t,c,k = interp.splrep(state.zmax_r,state.zmax, k = 2)
    f_bathy = interp.BSpline(t,c,k)
    z_int = f_bathy(r_int)

    plt.figure()
    plt.plot(state.zmax_r, state.zmax)
    plt.plot(r_int, z_int)
    plt.ylim(1200,0)

    plt.show()
    

