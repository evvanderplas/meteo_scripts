#! /usr/bin/env python

'''
sample script to demonstrate use of slice functionality
'''

import os,sys,glob
import numpy as np

levels,colors = None,None
gdir = '/nobackup_1/users/plas/temp'
gfile = os.path.join(gdir,'HARM_N25_201210250000_01200_GB')

print gfile, os.path.exists(gfile)

#grbs = pygrib.open(gfile)

import slice_class as gslice
import settings


par = 'soep'
if par == 'wind':
    u,v = gslice.gslice(),gslice.gslice()
    u.read_from_file(gfile,33)
    v.read_from_file(gfile,34)

    old = False # natuurlijk
    if old:
        myslice = u
        for l in myslice.data.keys():
            myslice.data[l] = np.sqrt(u.data[l]**2 + v.data[l]**2)
    else:
        myslice = (u**2 + v**2)**(1/2.)

elif par == 'soep':
    clwat,pcp,snow,graup = gslice.gslice(),gslice.gslice(),gslice.gslice(),gslice.gslice()
    clwat.read_from_file(gfile,76)
    pcp.read_from_file(gfile,181) #?
    snow.read_from_file(gfile,184) #?
    graup.read_from_file(gfile,201) #?

    ## calculate extinction 
    extinc = 144.7 * (1.2e3*clwat)**(0.88) \
        + 1.1 * (1.2e3*pcp)**(0.75) \
        + 10.4 * (1.2e3*snow)**(0.78) \
        + 2.4 * (1.2e3*graup)**(0.78) 
    myslice = 3912 * extinc**(-1) 
    levels = settings.VIS['plot_levels']
    colors = settings.VIS['colors']
    
elif par = 'TKE':
    myslice = gslice.gslice()
    myslice.read_from_file(gfile,200)

# set path
ll1 = (52.3,4.)
ll2 = (52.4,5.)
myslice.sline(ll1,ll2,40,verb=True)

# interpolate to the points in the path
myslice.intSliceVals()


# plot the slice
myslice.plot(levels = levels, colors = colors)
# plot a map to see where the slice is
myslice.plotmap(domain='nl')
