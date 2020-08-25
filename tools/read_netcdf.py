#! /usr/bin/env python

import os,sys,glob

import datetime as dt
import numpy as np
import netCDF4

ncdir      = '/net/bhw379/nobackup_1/users/plas/ncdat/'
ncfilename = 'KBYX_DTA_20120826_205100.nc'
ncfile = os.path.join(ncdir,ncfilename)
print ncfile, os.path.exists(ncfile)

ncobj = netCDF4.Dataset(ncfile,'r')
print 'In this file we have: \n',10*'='
for k in ncobj.variables.keys(): print k, ncobj.variables[k]

data  = ncobj.variables['dsa'][:]
lats  = ncobj.variables['lat'][:]
lons  = ncobj.variables['lon'][:]
datatime  = ncobj.variables['time'][:]

print 'Shape data,lat,lon:',data.shape,lats.shape,lons.shape
##contourf(lons,lats,data[0,:,:])

#nlats,nlons = np.meshgrid(lats,lons)
nlons,nlats = np.meshgrid(lons,lats)

ncobj.close()

## now make a plotinfo geofield instance:
import plotinfo
import pygrib

# make a radar "object":
rad = plotinfo.plotinfo(model='keysradar')
rad.values  = data[0,:,:]
rad.latlons = (nlats,nlons)
radardate = dt.datetime.fromtimestamp(datatime[0]) # convert from POSIX (seconds since 1-1-1970)
rad.date = radardate.strftime('%Y%m%d')
rad.time = radardate.strftime('%H%M')

# make a grib/forecast "object": (you can take any forecast file you want, offcourse)
gribdir  = '/net/bhw379/nobackup/users/plas/ectrans/'
gribfile = 'sst.grib'
gf = os.path.join(gribdir,gribfile)

gr = plotinfo.plotinfo(model='harmonie')
gr.read_from_grib(gf,par=11,leveltype='sfc',level=0,TR=0)

# resample radar unto the forecast grid
rad.resample(gr)
## in ipython: contourf(rad.values,np.arange(0,5,0.1))
