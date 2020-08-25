#!/usr/bin/env python

import pygrib

grbs = pygrib.open('hir_prec.grb')

for key in grbs[1].keys(): 
    for d in ('time','dat'):
        if d in key.lower():
            print key, grbs[1][key]
    #print key, grbs[1][key]

lats,lons = grbs[1].latlons()

import matplotlib, matplotlib.pyplot
from mpl_toolkits.basemap import Basemap # map functionality
#from matplotlib import cm as CM

fig = matplotlib.pyplot.figure()
axm = fig.add_subplot(1,1,1)

if 1:
    res = 'h'
    mer_int = 1
    par_int = 1
    pwidth  = 400.e3
    pheight = 400.e3

    bmap = Basemap(ax=axm,
                  projection='lcc', # 'lcc' is Lambert Conformal Conic projection
                  resolution='h',
                  lat_0=52.1,lon_0=4.7, # De Bilt
                  lat_1=45.,lat_2=55.,
                  width=pwidth, height=pheight
                  )

    x,y = bmap(lons,lats)

    bmap.drawcoastlines()
    bmap.drawcountries()

    print dir(bmap)

    del(bmap)

    print bmap
