#! /usr/bin/env python

'''
Simple example script of how to use *bmap* fast map plotting class
'''

import os,sys
import pygrib
import matplotlib.pyplot
from matplotlib import cm as CM

import bmap

if __name__ == '__main__':

    examplefile = './ex2012062100+012grib'
    grbs = pygrib.open(examplefile)
    for grb in grbs:
        if grb.indicatorOfParameter == 61 and grb.level == 457:
            values    = grb.values
            lats,lons = grb.latlons()

    # domain: standard is 'nl', also 'eur' (Europe), 'rot' (Rotterdam) and 'ijs' (IJsselmeer) are available
    # can be easily tweaked in bmap.py
    domain = 'eur'
    modelname = 'Harmonie'
    parametername = 'Precipitation (liquid)'
    lsmask_colour = 'black' # no need to set this

    ## Instantiate a figure (create a canvas to draw on)
    global fig
    fig  = matplotlib.pyplot.figure()
    ax   = fig.add_subplot(1,1,1)

    # set title
    ax.set_title('Example of {model} on domain {dom}'.format(model = modelname,dom = domain))

    ## make a basemap instance, set location, size, draw meridians etc
    plot_map = bmap.myMap(domain = domain, modelname = modelname)
    plot_map.xymap(latlons = (lats,lons), domain = domain, modelname = modelname)
    plot_map.set_axes(ax)
    plot_map.dress_up_map(domain = domain,lsmask_colour = lsmask_colour)

    pcplevels = (0.1,0.3,1,3,10,30,100)

    # Use your own colours...
    pcpcolors = ('white','lightgrey','grey','LightCoral','red','black')
    cmap      = None

    # ... or use something like this:
    #pcpcolors = None
    #cmap      = CM.get_cmap('jet')
    #cmap      = CM.get_cmap('hot_r')

    exarg     = 'both' # extending colors beyond the given levels, could be 'both, 'neither', see matplotlib doc

    pcont = plot_map.bmap.contourf(plot_map.x,plot_map.y,values,
                                   pcplevels,cmap=cmap,colors = pcpcolors,extend=exarg,ax=ax)

    # colorbar:
    cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.2f') #format='%.1e') 
    cb.set_label(parametername,fontsize=9)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(9)
    #if grib_info['draw_line']: 
    #    cb.add_lines(plcont)

    plotfile = 'test_bmap.png'
    fig.savefig(plotfile, dpi = 100)
