#! /usr/bin/env python

import os,sys
import numpy as np
import datetime as dt

import matplotlib.pyplot as plt
import matplotlib.cm as CM

import pygrib

from os.path import join as pjoin

import plotinfo

if __name__ == '__main__':

    print 'Calculate PV'

    model = 'VETREFO'
    gdir = '/nobackup_1/users/plas/temp/2010/07/14/06'
    gformat = 'HARM_N25_{ymdh}00_{lt}00_GB'

    model = 'VETSTD'
    gdir = '/nobackup_1/users/plas/temp/VETSTD'
    gformat = 'fc{ymdh}+{lt}grib'
    gformat_alt = 'fc{ymdh}+{lt}grib_md'

    dtg = dt.datetime(2010,7,14,6,0)
    
    pref = 1000
    levels = [1000, 950, 925, 900, 850, 800, 700, 600, 500, 400, 300, 200]
    leveltype = 'pl'
    levels = [60,55,50,45,40,35,30,20]
    leveltype = 'ml'

    pdom = 'eur'

    for lt in range(10,12):
        gfile = pjoin(gdir,gformat.format(ymdh=dtg.strftime('%Y%m%d%H'),lt = str(lt).zfill(3)))
        print gfile, os.path.exists(gfile)

        tpot = {}
        pl   = {}
        
        #for i,p in enumerate(pressl):
        for i,l in enumerate(levels):
            u = plotinfo.plotinfo(model='REFO')
            v = plotinfo.plotinfo(model='REFO')
            t = plotinfo.plotinfo(model='REFO')
            p = plotinfo.plotinfo(model='REFO')

            u = u.read_from_grib(gfile,par=33,level=l,leveltype=leveltype,case='REFO')
            v = v.read_from_grib(gfile,par=34,level=l,leveltype=leveltype,case='REFO')
            t = t.read_from_grib(gfile,par=11,level=l,leveltype=leveltype,case='REFO')
            p = p.read_from_grib(gfile+'_md',par= 1,level=l,leveltype=leveltype,case='REFO')

            vort = u.copy()
            planet_vort = 1/(24*60*60) * np.sin(vort.latlons[0])
            vort.values =  -1* ((u.values[1:,1:] - u.values[:-1,1:]) - (v.values[1:,1:] - v.values[1:,:-1]))/2500.
            vort.values = planet_vort + np.pad(vort.values,((1,0),(1,0)),mode='constant')
            
            vort.colormap = 'RdBu_r'
            vort.name = 'Vorticity'
            vort.plot_levels   = np.arange(-2.e-3,2.01e-3,1.e-4)
            vort.plot(domain=pdom,out='vort_plot_{dom}_{pl}_{ymdh}+{lt}.png'.format(ymdh = dtg.strftime('%Y%m%d%H'),
                                                                        lt=str(lt).zfill(3),
                                                                              pl = l,
                                                                              dom=pdom))

            #print vort.values.shape
            #sys.exit(0)

            pl[l]   = p.values
            pref = pl[levels[0]]

            t.values = t.values * np.power(pref/p.values,0.286) # pot temp

            # make previous T_pot available
            tpot[l] = t.values

            t.colormap = 'jet'
            t.name = 'Potential T'
            t.plot_levels   = np.linspace(int(t.values.min()),int(t.values.max()),20) #np.arange(290,311,1)
            t.plot(domain=pdom,out='pt_plot_{dom}_{pl}_{ymdh}+{lt}.png'.format(ymdh = dtg.strftime('%Y%m%d%H'),
                                                                         lt=str(lt).zfill(3),
                                                                         pl = l,
                                                                         dom=pdom))

            pv = vort.copy()
            if i > 0:
                pv.values = vort.values * (tpot[levels[i]] - tpot[levels[i-1]])/(pl[levels[i]]-pl[levels[i-1]])
                pv.name = 'Potential Vorticity'

                # calculate min, max, and due to resampling there are some nodata values (65535) to be omitted:
                pvmax   = np.percentile(abs(pv.values[:700,:700]),99)

                pv.plot_levels = np.linspace(-pvmax,pvmax,20)
                pv.plot(domain=pdom,out='pv_plot_{dom}_{pl}_{ymdh}+{lt}.png'.format(ymdh = dtg.strftime('%Y%m%d%H'),
                                                                                    lt=str(lt).zfill(3),
                                                                                    pl = l,
                                                                                    dom=pdom))


            #sys.exit(0)
            #break

