#! /usr/bin/env python

'''
sample script to demonstrate use of slice functionality
'''

import os,sys,glob
import numpy as np


gdir = '/nobackup_1/users/plas/temp'
gfile = os.path.join(gdir,'fc2012082600+014grib')

print gfile, os.path.exists(gfile)

#grbs = pygrib.open(gfile)

import slice_class as gslice

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

# set path
ll1 = (22.5,-80.)
ll2 = (25. ,-79.5)
myslice.sline(ll1,ll2,40,verb=True)

# interpolate to the points in the path
myslice.intSliceVals()

# plot the slice
myslice.plot()
# plot a map to see where the slice is
myslice.plotmap(domain='cuba')
