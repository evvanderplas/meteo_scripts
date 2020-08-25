#! /usr/bin/env python

import os,sys
import netCDF4
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

#
# See the notebook at :
# http://nbviewer.ipython.org/github/Unidata/tds-python-workshop/blob/master/netcdf-by-coordinates.ipynb
#

class Kdtree_fast(object):
    def __init__(self, ncfile, latvarname, lonvarname):
        self.ncfile = ncfile
        self.latvar = self.ncfile.variables[latvarname]
        self.lonvar = self.ncfile.variables[lonvarname]        

        if len(self.latvar.shape) == 1: 
            self.latvar = np.tile(self.latvar[:],(self.lonvar.shape[0],1))
            self.latvar = self.latvar.T
            self.lonvar = np.tile(self.lonvar[:],(self.latvar.shape[0],1))
        # Read latitude and longitude from file into numpy arrays
        rad_factor = np.pi/180.0 # for trignometry, need angles in radians
        self.latvals = self.latvar[:] * rad_factor
        self.lonvals = self.lonvar[:] * rad_factor
        self.shape = self.latvals.shape
        clat,clon = np.cos(self.latvals),np.cos(self.lonvals)
        slat,slon = np.sin(self.latvals),np.sin(self.lonvals)
        clat_clon = clat*clon
        clat_slon = clat*slon
        triples = zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat))
        self.kdt = cKDTree(triples)

    def query(self,lat0,lon0):
        rad_factor = np.pi/180.0 
        lat0_rad = lat0 * rad_factor
        lon0_rad = lon0 * rad_factor
        clat0,clon0 = np.cos(lat0_rad),np.cos(lon0_rad)
        slat0,slon0 = np.sin(lat0_rad),np.sin(lon0_rad)
        dist_sq_min, minindex_1d = self.kdt.query([clat0*clon0,clat0*slon0,slat0])
        #print np.unravel_index(minindex_1d, self.shape)
        iy_min, ix_min = np.unravel_index(minindex_1d, self.shape)

        return iy_min,ix_min

synops = open('allsynop.list','r')
latlons = []
for s in synops:
    try:
        sid,slat,slon = s.split()[0:3]
    except:
        print 'Line corrupt?',s
    latlons.append((float(slat),float(slon)))
#sys.exit(0)


testfile = '/nobackup_1/users/plas/ncdat/BULL/harm_2014043003.nc'

nc = netCDF4.Dataset(testfile,'r')
ns = Kdtree_fast(nc,'lat','lon')
ll = (52.0, 4.10)

dummy = []
for ll in latlons[1135:1200]:
    print 'Query {0}, {1}'.format(ll[0],ll[1])
    iy,ix = ns.query(ll[0],ll[1])
    #print 'Closest lat lon at :', iy,ix,ns.latvar[iy,ix], ns.lonvar[iy,ix]
    #print 'Then alatf there is ',nc.variables['alatf'][:,iy,ix]
    dummy.append((iy,ix)) #nc.variables['alatf'][:,iy,ix])

print 'dummy',len(dummy)

t2m = nc.variables['T2m']
print t2m[1,[d[0] for d in dummy],[d[1] for d in dummy]]
#nc.close()
