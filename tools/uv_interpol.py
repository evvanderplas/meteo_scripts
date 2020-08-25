#! usr/bin/env python

import os,sys,numpy
import pygrib

import interpol


# this gets executed when the script is run, not if it is imported:
if __name__ == '__main__':

    ## some example grib file:
    gfile = '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL/ex2013062400+003grib_sa'

    # open file
    grbs = pygrib.open(gfile)
    for grb in grbs:
        # parse grib messages, look for u10, v10
        if grb.indicatorOfParameter == 33 and grb.levelType == 'sfc' and grb.level == 10:
            u10       = grb.values
            lats,lons = grb.latlons()

        elif grb.indicatorOfParameter == 34 and grb.levelType == 'sfc' and grb.level == 10:
            v10 = grb.values

    # calculate, wind, wind direction
    wind = numpy.sqrt(u10**2 + v10**2)
    wdir = numpy.arctan2(v10,u10) /numpy.pi * 180

    #print wind,lats; sys.exit(1)

    # now take a single lat,lon to interpolate to:
    lattest,lontest = (52.375193,4.878892)
    intpffval = interpol.interp(wind,(lats,lons),(lattest,lontest),
                                method="bilinear")
    intpddval = interpol.interp(wdir,(lats,lons),(lattest,lontest),
                                method="bilinear")

    print intpffval,intpddval

    # take multiple lats,lons from a ASCII file :
    f = open('points.txt','r')
    flat,flon = [],[]
    for line in f:
        flat.append(float(line.split()[0]))
        flon.append(float(line.split()[1]))
    f.close()

    print flat,flon,len(flat)

    # make a simple data structure (dictionary): per point get lat,lon, ff and dd:
    intpdict = {}
    for i,lat,lon in zip(range(len(flat)),flat,flon):
        intpdict[i] = {}
        intpdict[i]['lat'] = lat
        intpdict[i]['lon'] = lon
        intpdict[i]['ff']  = interpol.interp(wind,(lats,lons),(lat,lon),
                                             method="bilinear")
        intpdict[i]['dd']  = interpol.interp(wdir,(lats,lons),(lat,lon),
                                             method="bilinear")
       
    print 'nr, lat,lon, wind and direction in degrees:'
    for p in intpdict.keys():
        print p,intpdict[p]['lat'],intpdict[p]['lon'],intpdict[p]['ff'],intpdict[p]['dd']
        
