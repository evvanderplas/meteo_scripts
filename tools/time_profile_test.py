#! /usr/bin/env python

import os,sys,pickle
import pygrib
import datetime as dt
import numpy as np
import matplotlib, matplotlib.pyplot
import matplotlib.cm as cm

from suf import suf2

import profile_new as profile
from profile_new import par_T, par_hum, par_u, par_v, par_vv, \
    par_TKE, par_PCP, par_GRP, par_SNW, par_CLI, par_CLW

'''
Usage: set a directory and a gribfile format, eg
    data_dir = '/nobackup/users/plas/temp'
    gribfilename = 'HARM_N25_{d}{st}00_0{lt}00_GB' (entering date, starttime and 2 digit leadtime)

This creates a 2d array of a number of parameters (in ml_pars,ml_par_names), and 
pickles it for future reference in 'profpickle.pkl'. 
When this file is not yet there, or the argument 'new' is given in the call of the script, 
it parses all the files, and dumps data in the file. Else it will just read the pickle file present. 

'''


## follow recipe in profile_new.py

def read_profile(gribname,date,st,platlon,params,paramnames,maxlead=24,nlevels = 60):

    ## needs a 'par_T' and a 'par_hum' to calculate z!!

    ## allocate arrays:
    ml_profile_values = {}
    for ml_par_name in paramnames: #ml_par_names:
        ml_profile_values[ml_par_name] = np.zeros((maxlead,nlevels))
    alt = np.zeros((maxlead,nlevels))

    # to make sure that interpolation is done only once:
    intpdone = False

    leadtimes = range(maxlead) # for now every hour default
    for st in (0,):
        for lt in leadtimes:
            gf = gribfilename.format(d=dtg.strftime('%Y%m%d'),st=suf2(st),lt=suf2(lt))
            gribfile = os.path.join(data_dir,gf); print gribfile, os.path.exists(gribfile)
            if not os.path.exists(gribfile): continue
            
            if not intpdone:
                # do a dummy read just to retrieve the lats,lons and a,b arrays
                print 'preparing interpolation'
                dummyvals,lats,lons = profile.get_grib_values(gribfile,params[0],
                                                              leveltype='sfc',level=0,TR=0)
            
                # do this just once for a given profile
                lli = profile.LatLonInterpolator(lats, lons)
                intpdone = True

            print lli,dir(lli)

            grbs = pygrib.open(gribfile)
            for grb in grbs:
                if (grb.levelType == 'ml'):
                    # TODO1: get these lat,lon values from the radiosonde data
                    # TODO2: use 2 grib files for time interpolation
                    lat_ml,lon_ml = platlon #52.3,4.2
                    for ml_par, ml_par_name in zip(params,paramnames): #ml_par_names):
                        if (grb.indicatorOfParameter == ml_par):
                            print 'leadtime {}, filling level {} for par {}'.\
                                format(lt,grb.level, ml_par_name)
                            ml_profile_values[ml_par_name][lt,grb.level-1] = \
                                lli.interpolate(grb.values,lat_ml,lon_ml)

            t_fl   = ml_profile_values['par_T'][lt,:]
            hum_fl = ml_profile_values['par_hum'][lt,:]
            pa_prov = profile.press_alt_provider(gribfile,lat_ml,lon_ml)
            alt[lt,:] = pa_prov.z_mid

            verb = True
            verb = False
            if verb:
                for i in range(pa_prov.nlevels+1):
                    if i<pa_prov.nlevels:
                        print "i=",i," p_hl = ",pa_prov.p_half[i],\
                            "phi_hl[i]=",pa_prov.phi_hl[i],\
                            " t_fl = ",pa_prov.T[i]," hum_fl = ",pa_prov.hum[i],\
                            " z_hl[i] = ",pa_prov.z_hl[i]
 
    return ml_profile_values,alt

## already in profile_new!!

#par_T     =  11 # temperature          (index   5 in the list) [K]
#par_phi_s =   6 # surface geopotential (index   9 in the list) [m^2/s^2]
#par_P_s   =   1 # surface pressure     (index  20 in the list) [pa]
#par_hum   =  51 # specific humidity    (index 160 in the list) [kg/kg]
# par_RH_s  =  52 # surface relative humidity (index 161/162 in the list)
#par_u     =  33 # u wind component     (index 163 in the list) [m/s]
#par_v     =  34 # v wind component     (index 164 in the list) [m/s]
#par_vv    =  40 # vert. velocity       (index 168 in the list) [not given]
par_TKE   =  200 # turbulent kinetic energy 

data_dir = '/nobackup/users/plas/temp'
modelname    = 'Harmonie'
gribfilename = 'HARM_N25_{d}{st}00_0{lt}00_GB'


dtg = dt.datetime(2012,11,20,0)
dtg = dt.datetime(2013,10,28,0) # storm
st = 0
leadtimes = range(0,49,1)
maxlead = 49
leadtimes = range(1,maxlead,1) # test
latlon = (52.3,4.2) # Schiphol
latlon = (51.964915,4.897662) # Cabauw
# make arrays:

ml_pars           = [ par_T,   par_hum,   par_u,   par_v,   par_vv]
ml_pars           = [ par_T,   par_hum,   par_u,   par_v,   par_vv, 
                      par_TKE,par_PCP,par_GRP,par_CLW,par_CLI]

ml_par_names      = ['par_T', 'par_hum', 'par_u', 'par_v', 'par_vv']
ml_par_names      = ['par_T', 'par_hum', 'par_u', 'par_v', 'par_vv',
                     'par_TKE','par_PCP','par_GRP','par_CLW','par_CLI']

nlevels = 60
levels = range(1,nlevels+1)


profpickle = 'prof.pkl'
if 'new' in sys.argv: print 'New!'; 
print 'Pickle exists: ',os.path.exists(profpickle)
#sys.exit(0)

if os.path.exists(profpickle) and 'new' not in sys.argv:
    f = open(profpickle,'r')
    ml_profile_values,alt = pickle.load(f)
    f.close()

else:
    print 'Parsing'; #sys.exit(1)
    gribfiles = os.path.join(data_dir,gribfilename)
    ml_profile_values,alt = read_profile(gribfiles,dtg,0,latlon,
                                         ml_pars,ml_par_names,maxlead=maxlead,nlevels = nlevels)

    f = open(profpickle,'wb')
    pickle.dump((ml_profile_values,alt),f)
    f.close()


if 0:
    make_plot = True
    if make_plot:

        fig = matplotlib.pyplot.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        new_fontsize = 6

        #            ax    title             xlabel              par
        plotdefs = [(ax1, 'profile of T',   'T [K]',            'par_T'),
                    (ax2, 'profile of hum', 'humidity [kg/kg]', 'par_hum'),
                    (ax3, 'profile of u',   'u [m/s]',          'par_u'),
                    (ax4, 'profile of vv',  'vertical velocity [m/s]', 'par_vv')]


        for ax, title, xlabel, par in plotdefs:
            ax.plot(ml_profile_values[par][lt,:],pa_prov.z_mid*1.e-3,marker='o', markersize=2)
            ax.set_title(title,fontsize=new_fontsize)
            ax.set_xlabel(xlabel,fontsize=new_fontsize)
            ax.set_ylabel('Altitude [km]',fontsize=new_fontsize)
            for l in ax.get_xticklabels()+ax.get_yticklabels():
                l.set_fontsize(fontsize=new_fontsize)
        
        plotfile = 'testprofile_{st}_{lt}.png'.format(st=st,lt=lt)
        fig.savefig(plotfile, dpi = 150)

else:

    sttime = dtg + dt.timedelta(hours=st)
    times  = [sttime + dt.timedelta(hours  = lt) for lt in leadtimes]

    figc = matplotlib.pyplot.figure()
    axc  = figc.add_subplot(111)
    #axc.contourf(leadtimes,ml_profile_values['par_T'])
    #axc.contourf(range(1,maxlead),alt[1,::-1].transpose(),ml_profile_values['par_vv'][1:,::-1].transpose(),np.arange(-.1,.2,.01))

    par = 'wind'
    par = 'par_CLW'

    if par == 'wind':
        u = ml_profile_values['par_u']
        v = ml_profile_values['par_v']
        ml_profile_values['wind'] = np.sqrt(u**2 + v**2)
        relevantmin = 0
        step = 1
    elif par == 'par_CLW':
        step = 3.e-5

    relevantmin = 0

    relevantmax = max(ml_profile_values[par][:,10:].ravel())
    relevantmax = max(ml_profile_values[par][:,:].ravel())
    cmap = cm.get_cmap(name='Reds')#, lut=None)

    pcont = axc.contourf(times,alt[1,::-1].transpose(),ml_profile_values[par][1:,::-1].transpose(),
                         np.arange(relevantmin,relevantmax,step))
    #pcont2 = axc.contourf(times,alt[1,::-1].transpose(),ml_profile_values['par_CLI'][1:,::-1].transpose(),
    #                      np.arange(relevantmin,relevantmax,step),cmap = cmap,alpha = 0.7)

    axc.set_xlabel('Time')

    axc.set_ylim(0,5.e3)
    axc.set_ylabel('Altitude (m)')
    axc.set_title('Profile of {m},{p} at {lat},{lon}, date {d}'\
                      .format(m=modelname,p=par,
                              lat = latlon[0],lon=latlon[1],
                              d=dtg.strftime('%Y%m%d')))

    # create colorbar
    cb = figc.colorbar(pcont,shrink=0.9, extend='both',format='%.2e') #format='%.1e') 
    cb.set_label(par,fontsize=9)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(9)

    # format dates
    figc.autofmt_xdate()

    # save
    plotfile = 'testprofile_c_{st}_{p}.png'.format(st=st,p=par)
    figc.savefig(plotfile, dpi = 150)
            
