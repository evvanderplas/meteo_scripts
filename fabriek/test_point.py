#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle
import time
import datetime as dt
import numpy
import matplotlib, matplotlib.pyplot
from mpl_toolkits.basemap import Basemap # map functionality
from matplotlib import cm as CM

import pygrib

#import h5tools,gribtools,finddate

class post_field:

    def __init__(self,name,
                 paramId,leveltype,level,
                 datetime,leadtime,
                 data):

        self.name = name
        self.paramId   = paramId
        self.leveltype = leveltype
        self.level     = level
        self.datetime  = datetime
        self.leadtime  = int(leadtime)
        self.data = data

class post_data:

    def __init__(self):
        self.fieldlist = []

    def add(self,name,
            paramId,leveltype,level,
            datetime,leadtime,
            data):

        self.fieldlist.append(post_field(name,
                                         paramId,leveltype,level,
                                         datetime,leadtime,
                                         data) )
        # check if already there?
        for item in self.fieldlist:
            #if field.paramId not in 
            #self.fieldlist.append(post_field(name,data))
            pass

    def get(self,paramId,leveltype,level,datetime,leadtime):
        
        for item in self.fieldlist:
            if item.paramId == paramId and item.leveltype == leveltype \
                    and item.level == level and item.datetime == datetime \
                    and item.leadtime == leadtime:
                
                return item.data

datadirjan = '/net/bhw177/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF'
rexp_ECJAN = 'fc*gribec'

gribec_files = sorted( glob.glob( os.path.join(datadirjan,rexp_ECJAN)) )
#for file in gribec_files: print os.path.split(file)[1]

if __name__ == '__main__':
    t0 = time.clock()

    import pickle
    import tools.read_stations as read_stations
    import get_closest_gridpoint.get_closest_gridpoint  as getpoint

    tsdump = './testpoint.pkl'
    if os.path.isfile(tsdump):
        f = open(tsdump,'r')
        ts_obj = pickle.load(f)
        f.close()

        print ts_obj.stationdatadict['240'].obsdict[11].get_fcdata()
    else:
        print 'creating ',tsdump

        grbfile = gribec_files[-1]
        print grbfile
        grbs = pygrib.open(grbfile)
        lats,lons = grbs[1].latlons()
        grbs.close

        stationfile = 'tools/stations.txt'
        stations = read_stations.read_stat_ascii(stationfile)
        for o in stations:
            print o.NAME, o.LAT, o.LON, #type(o.LAT), type(o.LON),
            g = getpoint.grid_searcher(lats,lons,method=2)
            a1,b1           = g.get_closest_point(float(o.LAT),float(o.LON))
            a2,b2,c2,d2     = g.get_four_surrounding_points(float(o.LAT),float(o.LON))
            print a1,b1
            print a2,b2,c2,d2

        grbfiles = gribec_files[-10:]
        import plottypes,sources
        import tools.timeseries as ts

        source = sources.ECJAN
        stationfile = 'tools/stations.txt'
        stations    = read_stations.stationdict(stationfile,grbfiles[0], model = source)

        ts_obj = ts.stationlist()
        for station in stations:
            stat = stations[station]
            ts_obj.add_station(stat['name'],stat['latlon'],stat['stationid'])

    

        outdir = './'
        for grib_file in grbfiles:

            grbindx=pygrib.index(grib_file,
                                 'indicatorOfParameter',
                                 'indicatorOfTypeOfLevel',
                                 'level')
            
            plottypes.point_extract(grbindx, source, stations, ts_obj, outdir = outdir)

        f = open(tsdump,'wb')
        pickle.dump(ts_obj,f)
        f.close()

        print ts_obj.stationdatadict['240'].obsdict[11].get_fcdata()
        #for stat in stations:
        #    print stations[stat]['name'],stations[stat]['latlon'],stations[stat][]
