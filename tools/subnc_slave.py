#! /usr/bin/env python

import os,sys
import pygrib
import numpy as np
import cPickle as pickle

gfile,startdate,lt = sys.argv[1],sys.argv[2],sys.argv[3]
print gfile,os.path.exists(gfile)

stationlist = pickle.load(open('statlist.pkl','rb'))

y,x = 150,122

ltrans = {'ml':109,
          'sfc':105,
          'pl':100,
          '103':103,
          }

valdict = {}

grbs = pygrib.open(gfile)
nrmessages = grbs.messages
infocom = np.zeros((nrmessages+1,5+max(stationlist.keys())),dtype=float)
for i,grb in enumerate(grbs):
    #print gfile[-7:-5],grb.indicatorOfParameter,grb.level
    ind,leveltype,level,tr,tab = grb.indicatorOfParameter,grb.levelType,grb.level,grb.timeRangeIndicator,grb.table2Version

    try:
        infocom[i,0] = ind
        infocom[i,1] = ltrans[leveltype]
        infocom[i,2] = level
        infocom[i,3] = tr
        infocom[i,4] = tab
        for s in stationlist:
            l,m = stationlist[s]['ind']
            #print l,m,grb.values.shape,s
            #print ind,leveltype,level,lt,s,grb.values[l,m]
            infocom[i,4+s] = grb.values[l,m]
    except:
        print "problem with ",ind,leveltype,level,tr,sys.exc_info()[0]
        sys.exit(0)

    if 0:
        if valdict.has_key(ind):
            if valdict[ind].has_key(leveltype):
                if valdict[ind][leveltype].has_key(level):
                    if valdict[ind][leveltype][level].has_key(tr):
                        pass
                    else:
                        valdict[ind][leveltype][level][tr] = grb.values[y,x]
grbs.close()


with open('profarr_{dtg}_{lt}.pkl'.format(dtg=startdate,lt=lt),"wb") as handle:
    pickle.dump(infocom,handle)
