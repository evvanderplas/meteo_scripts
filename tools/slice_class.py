#! /usr/bin/env python

import os,sys,pickle
import numpy as np
from scipy.interpolate import griddata
#import matplotlib
import pygrib

home = os.getenv('HOME')
sys.path.append(os.path.join(home,'python/tools'))
import interpol,bmap

def toLatlon(latlons,index):

    return [(latlons[0][i],latlons[1][i]) for i in index]

def do_kdtree(combined_x_y_arrays,points):

    '''
    Use a KDTree method to find the closest gridpoint (from stackoverflow)
    '''

    import scipy.spatial
    
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)

    return indexes

def find_four_neighbours(latlons,(xlat,xlon),verb=True):

    '''
    Find the indices of the four surrounding gridpoints to do the interpolation using the values at these corners
    Uses the lats,lons and a lat,lon of a point of interest
    '''

    # find nearest gridpoint:
    flats,flons = latlons[0].ravel(),latlons[1].ravel()
    dim0,dim1   = latlons[0].shape

    llarray = np.dstack([flats,flons])[0]
    res2    = do_kdtree(llarray,(xlat,xlon))
    i0,j0 = res2%dim0, (res2-res2%dim0)/dim0
    if verb: print 'KDTree: ',i0,j0,flats[res2],flons[res2]

    # determine if point is inside path: (uses pnpoly, presumably)
    from matplotlib.path import Path
    
    jn0,in0 = res2%dim0, (res2-res2%dim0)/dim0

    # define 4 paths on all sides of nearest neighbour:
    sq = {}
    sq['ll']   = [(in0,jn0),    (in0+1,jn0),  (in0+1,jn0+1),(in0,jn0+1)]
    sq['lr']   = [(in0-1,jn0),  (in0,jn0),    (in0,jn0+1),  (in0-1,jn0+1)]
    sq['ur']   = [(in0-1,jn0-1),(in0,jn0-1),  (in0,jn0),    (in0-1,jn0)]
    sq['ul']   = [(in0,jn0-1),  (in0+1,jn0-1),(in0+1,jn0),  (in0,jn0)]

    # loop over 4 adjacent grid cells:
    for corner in sq:
        path = sq[corner]
        llpath = Path(toLatlon(latlons,path))
        if verb: print 'Corner: ',corner,path,'\n',llpath

        # if the path contains the point of interest, return the indices:
        if llpath.contains_point((xlat,xlon)):
            if verb: print 'Succes! ',corner,path,llpath
            return path

    print 10*'*','No path found!',xlat,xlon

    # old and incorrect method:
    if flats[res2] < line_minlat:
        if flons[res2] < line_minlon:
            if verb: print 'linksonder'
            i01,j01 = i0,j0
            i02,j02 = i0+1,j0
            i03,j03 = i0+1,j0+1
            i04,j04 = i0,j0+1
        else:
            if verb: print 'rechtsonder'
            i01,j01 = i0-1,j0
            i02,j02 = i0,j0
            i03,j03 = i0,j0+1
            i04,j04 = i0-1,j0+1
    else:
        if flons[res2] < line_minlon:
            if verb: print 'linksboven'
            i01,j01 = i0-1,j0
            i02,j02 = i0,j0
            i03,j03 = i0,j0+1
            i04,j04 = i0-1,j0+1
        else:
            if verb: print 'rechtsboven'
            i01,j01 = i0-1,j0-1
            i02,j02 = i0,j0-1
            i03,j03 = i0,j0
            i04,j04 = i0-1,j0

    altindices = [(j01,i01),(j02,i02),(j03,i03),(j04,i04)]
    if verb: print 'Found: ',altindices

    return altindices

class gslice():

    def __init__(self,case = 'example'):
        self.case = case

    def copy(self):
        # make a (deep) copy
        cp = gslice()
        for attr in self.__dict__.keys():
            cp.__dict__[attr] = self.__dict__[attr]
        return cp

    def __add__(self,other):
        out = self.copy()
        #if hasattr(other,'__module__'):
        try:
            for l in out.data.keys():
                out.data[l] = self.data[l] + other.data[l]
        except:
            ## maybe int or float
            try:
                for l in out.data.keys():
                    out.data[l] = self.data[l] + other
            except:
                print "Ce n'est pas possible: +"; sys.exit(0)
        finally:
            return out

    def __sub__(self,other):
        out = self.copy()
        #if hasattr(other,'__module__'):
        try:
            for l in out.data.keys():
                out.data[l] = self.data[l] - other.data[l]
        except:
            ## maybe int or float
            try:
                for l in out.data.keys():
                    out.data[l] = self.data[l] - other
            except:
                print "Ce n'est pas possible: -"; sys.exit(0)
        finally:
            return out

    def __mul__(self,other):
        out = self.copy()
        #if hasattr(other,'__module__'):
        try:
            for l in out.data.keys():
                out.data[l] = self.data[l] * other.data[l]
        except:
            ## maybe int or float
            try:
                for l in out.data.keys():
                    out.data[l] = self.data[l] * other
            except:
                print "Ce n'est pas possible: *"; sys.exit(0)
        finally:
            return out

    def __pow__(self,other):
        out = self.copy()
        #if hasattr(other,'__module__'):
        try:
            for l in out.data.keys():
                out.data[l] = np.power(self.data[l],other)
        except:
            print "Ce n'est pas possible: pow ",other; sys.exit(0)
        finally:
            return out


    def read_from_file(self,gribfile,param,leveltype='ml'):

        '''
        Reading data from a gribfile
        '''

        print 'Reading ',gribfile
        try:
            f = open('slicedata.pkl','w')
        except:
            print 'reading new'

        self.file = gribfile

        import pygrib
        grbs = pygrib.open(gribfile)

        # adding some basics
        self.latlons = grbs[1].latlons()
        self.date    = grbs[1].dataDate
        self.time    = grbs[1].dataTime
        
        self.nrlevels = grbs[1].numberOfVerticalCoordinateValues/2
        pv = grbs[1].pv
        self.ab     = (pv[:self.nrlevels],pv[self.nrlevels:])
        #print self.ab

        self.param = param # might come in handy

        # read all levels 
        self.data = {}
        for grb in grbs:
            #print grb.indicatorOfParameter, grb.indicatorOfTypeOfLevel
            if grb.indicatorOfParameter == param and grb.indicatorOfTypeOfLevel == leveltype: #'hybrid':
                self.pname = grb.name
                l = grb.level; #print l,self.ab[0][l],self.ab[1][l]
                self.data[l] = grb.values
            #sys.exit(0)
        grbs.close()

    def sline(self,ll1,ll2,steps,verb=False):

        '''
        define the lats and lons of the line 
        (naive: |endpoint - startpoint|/nrofsteps, not a geodete), 
        store in class: self.line
        '''
        
        lat1,lon1 = ll1
        lat2,lon2 = ll2
        dlat,dlon =  (lat2-lat1)/steps, (lon2-lon1)/steps
        self.line = ([lat1 + s * dlat for s in range(steps+1)],
                     [lon1 + s * dlon for s in range(steps+1)])
        self.steps = steps

        if verb:
            print 'Line defined on'
            for lat,lon in zip(self.line[0],self.line[1]):
                print lat,lon

    def minimal_rect_around_line(self,verb=False):

        '''
        Draw a minimal rectangle around the points to minimise computing time for SciPy griddata
        Not too minimal, apparently, hence the +/- 2 in definition for i1,i2,j1,j2
        '''

        # Create a minimal rectangle of points around the line:
        #line_minlat,line_maxlat = min(self.line[0]), max(self.line[0])
        #line_minlon,line_maxlon = min(self.line[1]), max(self.line[1])

        # needs get_closest_gridpoint: hopefully to be replaced by some scipy routine someday
        #indicesmin = min(interpol.get_closest_index(self.latlons,(line_minlat,line_minlon)))
        #indicesmax = max(interpol.get_closest_index(self.latlons,(line_maxlat,line_maxlon)))

        #i1,i2 = max(0,indicesmin[0]-2),max(0,indicesmin[1]-2)
        #j1,j2 = min(dim0,indicesmax[0]+2),min(dim1,indicesmax[1]+2)

        altmin = find_four_neighbours(self.latlons,(self.line[0][0],self.line[1][0])) 
        altmax = find_four_neighbours(self.latlons,(self.line[0][-1],self.line[1][-1])) 
        altmin.extend(altmax)

        # avoid keyError:
        dim0,dim1 = self.latlons[0].shape

        i1,i2 = max(0,min([p[0] for p in altmin])), max(0,min([p[1] for p in altmin]))
        j1,j2 = min(dim0,max([p[0] for p in altmin])), min(dim1,max([p[1] for p in altmin]))

        # avoid keyError:
        dim0,dim1 = self.latlons[0].shape

        if verb:
            print indicesmin,indicesmax
            print 'i1,i2,j1,j2',i1,i2,j1,j2
            print 'minlat: ',i1,i2,self.latlons[0][i1,i2]
            print 'minlon: ',i1,i2,self.latlons[1][i1,i2]
            print 'maxlat: ',j1,j2,self.latlons[0][j1,j2]
            print 'maxlon: ',j1,j2,self.latlons[1][j1,j2]

        self.minrect = [i1,i2,j1,j2]
        
        return [i1,i2,j1,j2]

    def intValues(self,field):#,points = self.line):

        '''
        interpolate using SciPy: griddata(), default at points in self.line
        Uses self.latlons and some 2D field
        '''
        
        # recall minimal rectangle:
        i1,i2,j1,j2 = self.minrect
        
        x = self.latlons[1][i1:j1,i2:j2].flatten()
        y = self.latlons[0][i1:j1,i2:j2].flatten()
        z = field[i1:j1,i2:j2].flatten()

        # use scipy's griddata to interpolate:
        intval = griddata((x, y), z, (self.line[1], self.line[0]), method='linear')

        return intval
        
    def intSliceVals(self,verb=False):

        nrDataLevels = len(self.data.keys())
        self.sdata = np.zeros((nrDataLevels,self.steps+1),dtype=float)

        # recall minimal rectangle:
        self.minimal_rect_around_line()
        i1,i2,j1,j2 = self.minrect
        #return 1

        # do the interpolation for every level:
        for l in self.data.keys():
            intval = self.intValues(self.data[l])
            self.sdata[l-1,:] = intval 
            if verb: print self.sdata[l-1,:]

        if verb: print 'slice array:',self.sdata #.shape

    def compute_plevels(self,verb=False):

        self.mslp = np.zeros((self.steps+1),dtype = float)

        grbs = pygrib.open(self.file)
        for grb in grbs:
            if grb.indicatorOfParameter == 1  and grb.indicatorOfTypeOfLevel == 'sfc' and grb.level == 0:
                p = grb.values
        grbs.close()

        pvals = self.intValues(p)

        nrDataLevels = len(self.data.keys())
        self.pdata = np.zeros((nrDataLevels,self.steps+1),dtype=float)

        for s,ppt in enumerate(pvals):
            for l in self.data.keys():
                self.pdata[l-1,s] = self.ab[0][l] + self.ab[1][l] * ppt

        if verb:
            print 'Slice for pressure:'
            print self.pdata[:,3]
                
    def plot(self,levels = None,colors=None,colormap=None,outdir = './'):

        '''
        plot the vertical slice
        '''

        if levels is not None or colors is not None:
            pass
            # Not sure if this poses problems:
            #if len(levels) != len(colors): print 'levels should match colors:',levels,colors; #sys.exit(1)

        if not self.__dict__.has_key('leadtime'):
            self.leadtime = 0

        import matplotlib.pyplot as plt
        import matplotlib.cm as CM

        cmap = CM.get_cmap(colormap)

        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1)
        
        if levels == None: levels = len(self.data.keys())

        # do the plotting:
        #cnt = ax.contourf(self.pdata[::-1,:],self.sdata[::-1,:])
        cnt = ax.contourf(self.sdata[::-1,:],levels,cmap = colormap,extend='both')
        if self.__dict__.has_key('lines'):
            ax.contour(self.sdata[::-1,:],self.lines,colors='k')
            try:
                ax.contour(self.sdata[::-1,:],self.lines,colors='k')
            except:
                print 'Something awry in line def:',self.lines

        # adding stuff:
        cb = fig.colorbar(cnt,shrink=0.9, extend='both',format='%.2e') 
        cb.set_label(self.pname,fontsize=9)
        ax.set_xlabel('Steps along line')
        ax.set_ylabel('Model level')
        ax.set_title('{nam}, slice along ({lat1},{lon1}) - ({lat2},{lon2})'.format(nam=self.pname,lat1 = self.line[0][0],lon1 = self.line[1][0], lat2 = self.line[0][-1],lon2 = self.line[1][-1]))

        # save to file
        outfig = os.path.join(outdir,'slice_'+self.case+'_'+str(self.param)+'.png')
        outfig = os.path.join(outdir,'slice_{case}_{par}_{date}+{lt}.png'.format(case=self.case,par=self.param,date=self.date,lt=str(self.leadtime).zfill(3)))
        fig.savefig(outfig)
        print 'Created ',outfig

    def plotmap(self,domain = 'nl',colors=None,colormap='jet',ml = 60,lsmask_colour = 'black',outdir = './'):

        '''
        Create a map with one of the levels and the line along which the slice is taken indicated
        '''

        import matplotlib.pyplot as plt
        import matplotlib.cm as CM

        if not self.__dict__.has_key('leadtime'):
            self.leadtime = 0

        # create a figure
        fig  = plt.figure()
        ax   = fig.add_subplot(1,1,1)

        modelname = 'Harmonie'
        parametername = 'Parameter '+str(self.param)
        cmap = CM.get_cmap(colormap)

        # create the map
        plot_map = bmap.myMap(domain = domain, modelname = modelname)
        plot_map.xymap(latlons = (self.latlons[0],self.latlons[1]), domain = domain, modelname = modelname)
        plot_map.set_axes(ax)
        plot_map.dress_up_map(domain = domain,lsmask_colour = lsmask_colour)
        pcont = plot_map.bmap.contourf(plot_map.x,plot_map.y,self.data[ml],# here takes the last, maybe default to level?
                                       20,cmap=cmap,extend='both',ax=ax)

        # colorbar:
        cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.2f') #format='%.1e') 
        cb.set_label(parametername,fontsize=9)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(9)

        ax.set_title('Parameter '+str(self.param)+' (level '+str(ml)+'), * start, p end', fontsize=9)
        #xa,ya   = plot_map.bmap(4.878867,52.375345) #NL148
        #print xa,ya

        # calculate where the points are on the map and draw them
        xs,ys = plot_map.bmap(self.line[1],self.line[0]) 
        plot_map.bmap.plot(xs,ys,'o',color='white')

        # start and endpoint
        plot_map.bmap.plot(xs[0],ys[0],'*',color='white',markersize=15)
        plot_map.bmap.plot(xs[-1],ys[-1],'p',color='white',markersize=15)
        
        # save to file
        plotfile = os.path.join(outdir,'slice_map_'+self.case+'_'+str(self.param)+'.png')
        plotfile = os.path.join(outdir,'slice_map_{case}_{par}_{date}+{lt}.png'.format(case=self.case,par=self.param,date=self.date,lt=str(self.leadtime).zfill(3)))
        fig.savefig(plotfile, dpi = 100)
        print 'Created ',plotfile

if __name__ == '__main__':

    home = os.getenv('HOME')
    
    datadir = os.path.join(home,'HARP/Harp_sample/data')
    datadir = '/nobackup/users/plas'
    #datadir = os.path.join(home,'trans')
    gribfile = os.path.join(datadir,'fc','harm','ex2012062100+012grib')
    #gribfile = '/home/plas/data/radarver/fc2010071400+017grib'
    gribfile = os.path.join(datadir,'fc2013061806+012grib_37h12')

    myslice = gslice()
    myslice.read_from_file(gribfile,11)
    
    if 1:
        ll1 = (52.3,5.8)
        ll2 = (52.5,3.8)
        myslice.sline(ll1,ll2,30,verb=True)
        myslice.intSliceVals() #values()
        myslice.compute_plevels()
        myslice.plot()
        myslice.plotmap()
