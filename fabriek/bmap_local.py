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
#from tools.memory_inspector import report_mem_usage

# constant definitions
# ...

def suf3(h):
    if int(h)<10: return '00'+str(h)
    elif int(h)<100: return '0'+str(h)
    else: return str(h)

class myMap:

    nlmap  = {'alias':'nl',
              'file':'./nlmap.pkl',
              'res':'h',
              'latlon0':(52.1,4.7),
              'dim':(5.e5,5.e5)}
    eurmap = {'alias':'eur',
              'file':'./eurmap.pkl',
              'res':'i',
              'latlon0':(51.,3.),
              'dim':(20.e5,20.e5)}
    rotmap = {'alias':'rot',
              'file':'./rotmap.pkl',
              'res':'h',
              'latlon0':(52.,4.5),
              'dim':(1.2e5,1.2e5)}
    ijsmap = {'alias':'ijs',
              'file':'./ijsmap.pkl',
              'res':'h',
              'latlon0':(52.7,5.3),
              'dim':(1.e5,1.e5)}

    def __init__(self, domain = 'nl', modelname = 'harmonie', 
                 res = 'h',latlon0=(52.1,4.7), dim = (5.e5,5.e5)):

        maps = (self.nlmap, self.eurmap, self.rotmap, self.ijsmap) 

        if not hasattr(self.__class__,'pmap'):
            self.__class__.pmap = {}


        if not self.__class__.pmap.has_key((domain, modelname)):
            
            #check if pickle dump file exists
            map_index = None
            for (i,m) in enumerate(maps):
                if m['alias'] == domain:
                    map_index = i
                    
            if map_index is not None:
                mapdef = maps[map_index]
                if mapdef['alias'] == domain:
                    
                    if os.path.isfile(mapdef['file']): #open dumpfile, read bmap
                        print 'read map',str(domain)
                        
                        f = open(mapdef['file'],'r')
                        self.bmap = pickle.load(f)
                        f.close()

                    else: # construct instance
                        self.bmap = Basemap(
                            projection='lcc', 
                            resolution=mapdef['res'],
                            lat_0=mapdef['latlon0'][0],lon_0=mapdef['latlon0'][1],
                            lat_1=45.,lat_2=55.,
                            width=mapdef['dim'][0], height=mapdef['dim'][1]
                            )
                        
                        print 'write map ',domain
                        f = open(mapdef['file'],'wb')
                        pickle.dump(self.bmap,f)
                        f.close()                

            else: # map definition not recognized
                print 'Domain not pre-defined: ',domain
                try:
                    self.bmap = Basemap(
                        projection='lcc', 
                        resolution=res,
                        lat_0=latlon0[0],lon_0=latlon0[1],
                        lat_1=45.,lat_2=55.,
                        width=dim[0], height=dim[1]
                        )
                    print 'map was constructed, not written to file'
                except: # not enough or not the right arguments
                    'Define or modify domain in myMap class'
                    sys.exit(1)

            # finally: construct bmap instance
            self.__class__.pmap[(domain, modelname)] = self.bmap

        else: # read from self.__class__ :
            self.bmap = self.__class__.pmap[(domain, modelname)]

 
        self.domain = domain
        self.ax = "None"

        #self.xymap(file)
        self.cl  = None
        self.cnt = None
        self.par = None
        self.mer = None

        #self.dress_up_map(domain = domain)

    def set_axes(self,ax):

        self.ax = ax

    def dress_up_map(self,domain = None,
                 par_bot = 45, par_top = 60, par_int = 1, 
                 mer_bot = -1, mer_top = 10, mer_int = 1,
                 labelfont = 6,lsmask_colour = 'black'):

        if domain in ('nl','rot','ijs'): 
            par_bot, par_top, par_int = 45,60, 1 
            mer_bot, mer_top, mer_int = -1,10, 1
        elif domain == 'eur': 
            par_bot, par_top, par_int =  40,60, 2 
            mer_bot, mer_top, mer_int = -16,21, 2
        else:
            par_bot, par_top, par_int = 45,60, 1 
            mer_bot, mer_top, mer_int = -1,10, 1
           
        try:
            for key in self.cl:  del(self.cl[key])
            for key in self.cnt: del(self.cnt[key])
            for key in self.par: del(self.par[key])
            for key in self.mer: del(self.mer[key])
        except:
            pass

        #print 'to draw coastlines etc',time.clock() - t0, int(time.clock() - t0)*'='
        #print 'drawing coastlines ', domain
        self.cl  = self.bmap.drawcoastlines(ax=self.ax,color = lsmask_colour)
        self.cnt = self.bmap.drawcountries(ax=self.ax,color = lsmask_colour)
        #print self.cl


        #print 'to draw parallels etc',time.clock() - t0, int(time.clock() - t0)*'='
        self.par = self.bmap.drawparallels(numpy.arange(par_bot,
                                                   par_top,
                                                   par_int),
                                      labels=[1,1,0,0])
        self.mer = self.bmap.drawmeridians(numpy.arange(mer_bot,
                                                   mer_top,
                                                   mer_int),
                                      labels=[0,0,0,1])

        #print 'new font size in labels:',time.clock() - t0, int(time.clock() - t0)*'='
        new_fontsize = labelfont
        for key in self.par.keys():
            line,text = self.par[key]
            for txt in text:
                txt.set_fontsize(new_fontsize)
        for key in self.mer.keys():
            line,text = self.mer[key]
            for txt in text:
                txt.set_fontsize(new_fontsize)
                
    def xymap(self, latlons = [], files = [], domain = 'nl', modelname = 'harmonie'):

        '''Mapping the lats and lons to actual locations on the map '''


        #print '==> xymap start: ',
        #report_mem_usage()
        
        # create dict of possible x,y
        if not hasattr(self.__class__,'map_xy'):
            self.__class__.map_xy = {}

        if self.__class__.map_xy.has_key((domain, modelname)):
            # read x,y from class memory
            (self.x, self.y) = self.__class__.map_xy[(domain, modelname)]

        elif not self.__class__.map_xy.has_key((domain, modelname)):
            try:
                lats,lons = latlons
                self.x,self.y = self.bmap(lons,lats)

                # write to class memory:
                self.__class__.map_xy[(domain, modelname)] = (self.x,self.y)

            except:
                try:
                    # open file, extract lats,lons and map to x,y
                    grbs = pygrib.open(files[0])
                    lats,lons = grbs[1].latlons()

                    self.x,self.y = self.bmap(lons,lats)

                    # write to class memory:
                    self.__class__.map_xy[(domain, modelname)] = (self.x,self.y)
                
                except:
                    print 'no files in specified location',domain,modelname,files
                    sys.exit(1)

        #print '==> xymap end:   ',
        #report_mem_usage()

    def testmap(self,fig,ax,var,param,level,
                plotlevels = range(10),draw_line = False, 
                cmap = None,colors = None,
                name = 'Variable',date = 20040514, leadtime = 0, outdir = './'):

        '''A plotting method, not generally used,
        needs a figure, axis instance and some plotting information
        '''

        levels = plotlevels
        if cmap == None and colors == None: 
            cmap = CM.get_cmap('jet')
        elif cmap is not None:
            cmap.set_over(color='k',alpha = None)
            colors = None
        else:
            pass
        
        leadtime = suf3(leadtime)

        ax.set_title(name+' (parameter '+str(param)+', level '+str(level)+') at '+str(date)+'+'+str(leadtime) )
        pcont = self.bmap.contourf(self.x,self.y,var,levels,cmap=cmap,colors = colors,ax=self.ax)
        if draw_line: 
            plcont = self.bmap.contour(self.x,self.y,var,levels,ax=self.ax)

        # colorbar:
        cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.1e') 
        cb.set_label(name,fontsize=9)

        # naming and saving
        plotfile = os.path.join(outdir,
                                'test_'+self.domain+\
                                    '_'+str(param)+'_'+str(level)+'_'\
                                    +str(date)+'+'+str(leadtime)+'.png')
        sf = fig.savefig(plotfile)
        print 'created plot: '+plotfile
        
       #clean up
        for coll in pcont.collections:
            coll.remove()
    
        try:
            for coll in plcont.collections:
                coll.remove()
        except:
            pass

        # delete colorbar
        cbcoll = fig.get_axes()[1]
        cbcoll.collections = []
        fig.delaxes(cbcoll)
        fig.subplots_adjust(right=0.90)

        # return fig #?
