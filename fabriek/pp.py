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
        self.leadtime  = leadtime
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

    def get_specific(self,paramId,leveltype,level,datetime,leadtime):

        for item in self.fieldlist:
            if item.paramId == paramId and item.leveltype == leveltype \
                    and item.level == level and item.datetime == datetime \
                    and item.leadtime == leadtime:
                
                return item.data

    def get_all(self,paramId,leveltype,level):

        datalist = []
        for item in self.fieldlist:
            if item.paramId == paramId \
                    and item.leveltype == leveltype \
                    and item.level == level:

                datalist.append( item.data )

        return datalist

