#! /usr/bin/env python

import os,sys,pygrib
import datetime as dt
import numpy as np

from suf import suf2

import profile_new as profile

## follow recipe 
par_T     =  11 # temperature          (index   5 in the list) [K]
par_phi_s =   6 # surface geopotential (index   9 in the list) [m^2/s^2]
par_P_s   =   1 # surface pressure     (index  20 in the list) [pa]
par_hum   =  51 # specific humidity    (index 160 in the list) [kg/kg]
# par_RH_s  =  52 # surface relative humidity (index 161/162 in the list)
par_u     =  33 # u wind component     (index 163 in the list) [m/s]
par_v     =  34 # v wind component     (index 164 in the list) [m/s]
par_vv    =  40 # vert. velocity       (index 168 in the list) [not given]
par_TKE   =  200 # turbulent kinetic energy 

data_dir = '/net/bens03/harmonie_data'
data_dir = '/net/bhw379/nobackup/users/plas/temp'
data_dir = '/net/bhw277/nobackup/users/westrhen/radiosonde/nwp'

# can be looped
dtg = dt.datetime(2012,11,20,0)
dtg = dt.datetime(2013,10,24,0)
st = '00'
lt = '00'
gribfilename = 'HARM_N25_{d}{st}00_0{lt}00_GB'

# lat,lon where profile should be taken:
lat_ml,lon_ml = 52.3,4.2
lat_ml,lon_ml = 52.1009,5.1776

# make arrays:
ml_pars           = [ par_T,   par_hum,   par_u,   par_v,   par_vv, par_TKE]
ml_pars           = [ par_T,   par_hum,   par_u,   par_v,   par_vv]
ml_par_names      = ['par_T', 'par_hum', 'par_u', 'par_v', 'par_vv','par_TKE']
ml_par_names      = ['par_T', 'par_hum', 'par_u', 'par_v', 'par_vv']

nlevels = 60
levels = range(1,nlevels+1)

# initialise
ml_profile_values = {}
for ml_par_name in ml_par_names:
    ml_profile_values[ml_par_name] = np.zeros(nlevels)

gf = gribfilename.format(d=dtg.strftime('%Y%m%d'),st=suf2(st),lt=suf2(lt))
gribfile = os.path.join(data_dir,gf); print gribfile, os.path.exists(gribfile)

# do a dummy read just to retrieve the lats,lons and a,b arrays
print 'preparing interpolation'
dummyvals,lats,lons = profile.get_grib_values(gribfile,par_T,
                                              leveltype='sfc',level=0,TR=0)

# do this just once for a given profile
lli = profile.LatLonInterpolator(lats, lons)

# extract info from gribfile
grbs = pygrib.open(gribfile)
for grb in grbs:
    if (grb.levelType == 'ml'):
        for ml_par, ml_par_name in zip(ml_pars,ml_par_names):
            if (grb.indicatorOfParameter == ml_par):
                #print 'filling level {} for par {}'.format(grb.level, ml_par_name)
                ml_profile_values[ml_par_name][grb.level-1] = lli.interpolate(grb.values,lat_ml,lon_ml)

# calculate height (z, m) model levels:
t_fl   = ml_profile_values['par_T']
hum_fl = ml_profile_values['par_hum']
pa_prov = profile.press_alt_provider(gribfile,lat_ml,lon_ml)

# now ml_profile_values[par] (eg par = par_T for temperature) contains the profile information
# and pa_prov.z_mid contains the height (of the mid-levels) in meters
# In ipython --pylab one could
import matplotlib.pyplot as plt
fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
new_fontsize = 6
for par in ml_par_names:
    # make plot
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    # plot profile
    ax.plot(ml_profile_values[par],pa_prov.z_mid*1.e-3)
    # add some features
    ax.set_xlabel(par,fontsize=new_fontsize)
    ax.set_ylabel('Altitude [km]',fontsize=new_fontsize)
    # save 
    figname = 'prof_{par}_{dtg}_{st}_{lt}.png'.format(par = par,dtg= dtg.strftime('%Y%m%d'),st=st,lt=lt)
    fig.savefig(figname,dpi=200)
