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

    if 1:

        param_list = [[61,105,457,[1,2,3,4]],[61,105,456,numpy.arange(0,5e-3,5e-4)],[11,105,2,numpy.arange(260,310,1)]]
        for l in range(len(param_list)): 
            if param_list[l][1] == 105  : 
                print l,param_list[l]
                param_list[l][1] = 'sfc'
            elif param_list[l][1] == 109: param_list[l][1] = 'ml'

        searchlist = [ param_list[item][0:3] for item in range(len(param_list))]
        levellist  = [ param_list[item][3]   for item in range(len(param_list))]
        print searchlist,levellist;#sys.exit(1)

        a = post_data()
        print a

        for file in gribec_files[0:3]:

            print file
            grbindx=pygrib.index(file,'indicatorOfParameter','indicatorOfTypeOfLevel','level')
            print grbindx

            #selected_grbs=grbindx.select(indicatorOfParameter=61,indicatorOfTypeOfLevel='sfc',level=[456,457]) #,typeOfLevel='isobaricInhPa',level=500)
            #for grb in selected_grbs: 
            #    print grb
            #    print grb.indicatorOfTypeOfLevel            


            grbs = pygrib.open(file)
            print grbs(indicatorOfParameter=61,indicatorOfTypeOfLevel='sfc',level=[456,457])
            sys.exit(1)
            for grb in grbs:
                #if (grb.indicatorOfParameter,grb.indicatorOfTypeOfLevel,grb.level) in param_list:
                #if grb.indicatorOfParameter in (11,61) and grb.level in (2,457):
                if [grb.indicatorOfParameter,grb.indicatorOfTypeOfLevel,grb.level] in searchlist:
                    print grb
                    print grb.name
                    print grb.indicatorOfTypeOfLevel
                    ind = searchlist.index([grb.indicatorOfParameter,grb.indicatorOfTypeOfLevel,grb.level])

                    w     = grb.values
                    param = grb.indicatorOfParameter
                    ltype = grb.indicatorOfTypeOfLevel
                    level = grb.level
                    datetime = grb.dataDate
                    dtime    = grb.dataTime
                    leadtime = grb.stepRange
                    name  = grb.name
                    #lats,lons = grb.latlons()
                    
                    print w.min(), w.max(),levellist[ind]
                          
                    
                    a.add(name,
                          param,ltype,level,datetime,leadtime,
                          w)
                    print 'added', name,param,ltype,level,datetime,leadtime

        #def get(self,paramId,leveltype,level,datetime,leadtime):
        print dir(a)
        for item in a.fieldlist:
            print item.paramId,item.leveltype,item.level, item.datetime, item.leadtime
            print type(item.paramId),type(item.leveltype),type(item.level), type(item.datetime), type(item.leadtime)
        #2 metre temperature 11 sfc 2 20120304 3
        #b = a.get(11,'sfc',2,20120304,3)
        b = a.get(11,'sfc',2,20120225,17)
        #print b, dir(b)
        print b.shape, b.max()
