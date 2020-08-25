#! /usr/bin/env python

import os,sys,pickle
import numpy as np
from   scipy.interpolate import Rbf,griddata



# import tool to find nearest gridpoints, calculate distances on a sphere etc
home = os.getenv('HOME')
toolpath = os.path.join(home,'python','tools'); print 'Using path', toolpath
sys.path.append(os.path.join(home,'python','tools'))
import get_closest_gridpoint  as getpoint

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
    i0,j0 = res2%dim0, (res2-res2%dim0)/dim0 # Harmonie?
    i0,j0 = res2%dim1, (res2-res2%dim1)/dim1 # D11, ECMWF (?)
    if verb: 
        print 'KDTree: ',i0,j0,flats[res2],flons[res2],latlons[0][j0,i0],latlons[1][j0,i0]
        #print 'goodindex?',np.where(latlons[0] == flats[res2])
        #print 'goodindex?',np.where(latlons[1] == flons[res2])

    # determine if point is inside path: (uses pnpoly, presumably)
    from matplotlib.path import Path
    
    jn0,in0 = res2%dim0, (res2-res2%dim0)/dim0
    jn0,in0 = res2%dim1, (res2-res2%dim1)/dim1

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
        if verb: 
            print 'For point ',xlat,xlon, 'from nn ',jn0,in0
            print 'Corner: ',corner,path,'\n',llpath

        # if the path contains the point of interest, return the indices:
        if llpath.contains_point((xlat,xlon)):
            if verb: 
                print 'Succes! ',corner,path,llpath
                #for ll in list(llpath):
                #    print 'So backsubstitute in lats,lons ',ll,latlons[0][ll],latlons[1][ll]
            return path

    print 10*'*','No path found!',xlat,xlon


class reinterpol():

    def __init__(self,gridll = None,latlonlist = []):

        self.stationll = latlonlist
        self.latlons   = gridll
        
        self.indices = {} # dict of indices in the gridded field
        # dicts of values at the corners (for each point of interest)
        self.vala    = {} 
        self.valb    = {}
        self.valc    = {}
        self.vald    = {}
        
        self.xa,self.xb,self.xc,self.xd = {},{},{},{}
        self.ya,self.yb,self.yc,self.yd = {},{},{},{}
        self.distab,self.distbc,self.distcd,self.distda = {},{},{},{}
        self.xcorr, self.ycorr = {},{}


        self.yab,self.ycd = {},{}

        self.data    = None
        
    def compute_dist(self,verb=True):
    
        '''
        get all distances in "square" around point of interest:

        d - c    a - b
        | t | or |   | 
        a - b    d - c
        '''
        lats,lons = self.latlons
        
        for s,ll in enumerate(self.stationll):
            latt,lont = ll
            if verb: print 'point of interest (lat,lon) ',s,latt,lont
            g        = getpoint.grid_searcher(lats,lons,method=2)
            a,b,c,d  = g.get_four_surrounding_points(latt,lont)
            index    = (a,b,c,d)
            print 'Index:',index
            self.indices[s] = ((index[0],index[2]),
                               (index[0],index[3]),
                               (index[1],index[3]),
                               (index[1],index[2]) )
            if verb: print 'Indices of 4 surrounding points: ',self.indices[s]

            lata,lona = lats[self.indices[s][0]],lons[self.indices[s][0]]
            latb,lonb = lats[self.indices[s][1]],lons[self.indices[s][1]]
            latc,lonc = lats[self.indices[s][2]],lons[self.indices[s][2]]
            latd,lond = lats[self.indices[s][3]],lons[self.indices[s][3]]

            self.xa[s]     = getpoint.get_distance_grcircle(lata,lona,lata,lont)
            self.xb[s]     = getpoint.get_distance_grcircle(latb,lonb,latb,lont)
            self.xc[s]     = getpoint.get_distance_grcircle(latc,lonc,latc,lont)
            self.xd[s]     = getpoint.get_distance_grcircle(latd,lond,latd,lont)
            
            self.ya[s]     = getpoint.get_distance_grcircle(lata,lona,latt,lona)
            self.yb[s]     = getpoint.get_distance_grcircle(latb,lonb,latt,lonb)
            self.yc[s]     = getpoint.get_distance_grcircle(latc,lonc,latt,lonc)
            self.yd[s]     = getpoint.get_distance_grcircle(latd,lond,latt,lond)

            self.distab[s] = getpoint.get_distance_grcircle(lata,lona,latb,lonb)
            self.distbc[s] = getpoint.get_distance_grcircle(latc,lonc,latb,lonb)
            self.distcd[s] = getpoint.get_distance_grcircle(latc,lonc,latd,lond)
            self.distda[s] = getpoint.get_distance_grcircle(lata,lona,latd,lond)

            #self.yab[s] = self.ya[s] * self.xb[s] / self.distab[s] + self.yb[s] * self.xa[s] / self.distab[s]  
            #self.ycd[s] = self.yd[s] * self.xc[s] / self.distab[s] + self.yc[s] * self.xd[s] / self.distcd[s] 

            self.xcorr[s] = (self.distab[s] / (self.xa[s] + self.xb[s]))
            self.ycorr[s] = (self.distda[s] / (self.ya[s] + self.yd[s]))

            if verb: 
                print 'check x (xa+xb,xc+xd)',self.xa[s] + self.xb[s],self.xc[s] + self.xd[s]
                print 'check y (ya+yd,yb+yc)',self.ya[s] + self.yd[s],self.yb[s] + self.yc[s]
                print 'check corr x',(self.xa[s] + self.xb[s])*(self.xcorr[s])
                print 'check corr y',(self.ya[s] + self.yd[s])*(self.ycorr[s])
                print 'check dist (ab,bc,cd,da)',self.distab[s],self.distbc[s],self.distcd[s],self.distda[s]
                
    def multi_intp(self,valuedict,verb=True):

        self.data = np.zeros([len(valuedict.keys()),len(self.stationll)])

        for s in range(len(self.stationll)):

            for key in valuedict.keys():
            
                a = valuedict[key][self.indices[s][0]]
                b = valuedict[key][self.indices[s][1]]
                c = valuedict[key][self.indices[s][2]]
                d = valuedict[key][self.indices[s][3]]
                #print self.indices[s]
        
                # bilinear:
                #pab = a * self.xb[s] / self.distab[s] + b * self.xa[s] / self.distab[s]  
                #pcd = d * self.xc[s] / self.distab[s] + c * self.xd[s] / self.distcd[s] 

                wa = self.yd[s] * self.xb[s] * self.xcorr[s] * self.ycorr[s]/ (self.distda[s] * self.distab[s])
                wb = self.yc[s] * self.xa[s] * self.xcorr[s] * self.ycorr[s]/ (self.distbc[s] * self.distab[s])
                wc = self.yb[s] * self.xd[s] * self.xcorr[s] * self.ycorr[s]/ (self.distbc[s] * self.distcd[s])
                wd = self.ya[s] * self.xc[s] * self.xcorr[s] * self.ycorr[s]/ (self.distda[s] * self.distcd[s])
                print 'weights:',wa,wb,wc,wd,',sum:',wa+wb+wc+wd
                
                pvalue = a * wa + b * wb + c * wc + d * wd
                if verb:
                    print  a,wa,'\n',b,wb,'\n',c,wc,'\n',d,wd,'\n',10*'-'
                    print pvalue,(wa+wb+wc+wd),'!=',(a+b+c+d)/4,(a+b+c+d)*(wa+wb+wc+wd)/4
                self.data[key-1,s] = pvalue

        return self.data 
# end class
         
def get_closest_index(latlons,stationll):

    """
    Returns the indices of the four surrounding gridpoints in _latlons_ around _stationll_
    """
    import get_closest_gridpoint  as getpoint

    lats,lons = latlons
    latt,lont = stationll
    print 'point of interest (lat,lon) ',latt,lont

    g        = getpoint.grid_searcher(lats,lons,method=2)
    a,b,c,d  = g.get_four_surrounding_points(latt,lont)
    index    = (a,b,c,d)
    print 'Index:',index
    indices = ((index[0],index[2]),
               (index[0],index[3]),
               (index[1],index[3]),
               (index[1],index[2]))
    print 'Indices of 4 surrounding points: ',indices

    return indices

def interp(values,latlons,stationll,indices=None,method = 'bilinear'):

    import get_closest_gridpoint as getpoint

    lats,lons = latlons
    latt,lont = stationll

    try:
        f        = open('line_intp.pkl')
        linedict = pickle.load(f)
    except:
        print 'no predefined coefficients available, compute'
        
    if indices is None:
        print 'point of interest (lat,lon) ',latt,lont
        g        = getpoint.grid_searcher(lats,lons,method=2)
        a,b,c,d  = g.get_four_surrounding_points(latt,lont)
        index    = (a,b,c,d)
        print 'Index:',index
        indices = ((index[0],index[2]),
                   (index[0],index[3]),
                   (index[1],index[3]),
                   (index[1],index[2]))
        print 'Indices of 4 surrounding points: ',indices

    elif indices is not None and len(indices) < 4:
        return values[indices[0],indices[1]]

    if 1:
        a = values[indices[0]]
        b = values[indices[1]]
        c = values[indices[2]]
        d = values[indices[3]]
        #print indices,a,b,c,d

        
        lata,lona = lats[indices[0]],lons[indices[0]]
        latb,lonb = lats[indices[1]],lons[indices[1]]
        latc,lonc = lats[indices[2]],lons[indices[2]]
        latd,lond = lats[indices[3]],lons[indices[3]]

        if method == 'bilinear':
            latt,lont = stationll
            
            xa     = getpoint.get_distance_grcircle(lata,lona,lata,lont)
            xb     = getpoint.get_distance_grcircle(latb,lonb,latb,lont)
            xc     = getpoint.get_distance_grcircle(latc,lonc,latc,lont)
            xd     = getpoint.get_distance_grcircle(latd,lond,latd,lont)

            ya     = getpoint.get_distance_grcircle(lata,lona,latt,lona)
            yb     = getpoint.get_distance_grcircle(latb,lonb,latt,lonb)
            yc     = getpoint.get_distance_grcircle(latc,lonc,latt,lonc)
            yd     = getpoint.get_distance_grcircle(latd,lond,latt,lond)
            #yab    = getpoint.get_distance_grcircle(lata,lona,latt,lona)

            distap =  getpoint.get_distance_grcircle(lata,lona,latt,lont)
            distbp =  getpoint.get_distance_grcircle(latb,lonb,latt,lont)
            distcp =  getpoint.get_distance_grcircle(latc,lonc,latt,lont)
            distdp =  getpoint.get_distance_grcircle(latd,lond,latt,lont)

            distab = getpoint.get_distance_grcircle(lata,lona,latb,lonb)
            distbc = getpoint.get_distance_grcircle(latc,lonc,latb,lonb)
            distcd = getpoint.get_distance_grcircle(latc,lonc,latd,lond)
            distda = getpoint.get_distance_grcircle(lata,lona,latd,lond)

            # bilinear:
            pab = a * xb / distab + b * xa / distab  
            pcd = d * xc / distab + c * xd / distcd 

            yab = ya * xb / distab + yb * xa / distab  
            ycd = yd * xc / distab + yc * xd / distcd 
            pvalue = pab * ycd / distda + pcd * yab / distda

            if 0: # debug
                print a,b,xb,xa, distab
                print lata,latb,ya,yb,yab
                print latc,latd,yc,yd,ycd
                print pab,pcd,pvalue

            return pvalue

        else: # return average (?)
            return (a+b+c+d)/4.

def multi_interp(listdictvalues,latlons,stationlllist,indices=None,method = 'bilinear'):

    ''' 
    assume a number of fields to be interpolated, and a number of points
    Gather the interpolation characteristics to reuse(!)
    '''
    
    import get_closest_gridpoint  as getpoint

    lats,lons = latlons
    #latt,lont = stationll

    try:
        f        = open('line_intp.pkl')
        linedict = pickle.load(f)
    except:
        print 'no predefined coefficients available, compute'
        
        #if indices is None:
        indices = {} 
        for s,ll in enumerate(stationlllist):
            latt,lont = ll
            print 'point of interest (lat,lon) ',latt,lont
            g        = getpoint.grid_searcher(lats,lons,method=2)
            a,b,c,d  = g.get_four_surrounding_points(latt,lont)
            index    = (a,b,c,d)
            print 'Index:',index
            indices[s] = ((index[0],index[2]),
                          (index[0],index[3]),
                          (index[1],index[3]),
                          (index[1],index[2]))
            print 'Indices of 4 surrounding points: ',indices[s]

            lata,lona = lats[indices[0]],lons[indices[s][0]]
            latb,lonb = lats[indices[1]],lons[indices[s][1]]
            latc,lonc = lats[indices[2]],lons[indices[s][2]]
            latd,lond = lats[indices[3]],lons[indices[s][3]]

            if method == 'bilinear':
                xa,xb,xc,xd = {},{},{},{}
                ya,yb,yc,yd = {},{},{},{}
                distab,distbc,distcd,distda = {},{},{},{}
            
                xa[s]     = getpoint.get_distance_grcircle(lata,lona,lata,lont)
                xb[s]     = getpoint.get_distance_grcircle(latb,lonb,latb,lont)
                xc[s]     = getpoint.get_distance_grcircle(latc,lonc,latc,lont)
                xd[s]     = getpoint.get_distance_grcircle(latd,lond,latd,lont)

                ya[s]     = getpoint.get_distance_grcircle(lata,lona,latt,lona)
                yb[s]     = getpoint.get_distance_grcircle(latb,lonb,latt,lonb)
                yc[s]     = getpoint.get_distance_grcircle(latc,lonc,latt,lonc)
                yd[s]     = getpoint.get_distance_grcircle(latd,lond,latt,lond)
                #yab[s]    = getpoint.get_distance_grcircle(lata,lona,latt,lona)

                #distap =  getpoint.get_distance_grcircle(lata,lona,latt,lont)
                #distbp =  getpoint.get_distance_grcircle(latb,lonb,latt,lont)
                #distcp =  getpoint.get_distance_grcircle(latc,lonc,latt,lont)
                #distdp =  getpoint.get_distance_grcircle(latd,lond,latt,lont)

                distab[s] = getpoint.get_distance_grcircle(lata,lona,latb,lonb)
                distbc[s] = getpoint.get_distance_grcircle(latc,lonc,latb,lonb)
                distcd[s] = getpoint.get_distance_grcircle(latc,lonc,latd,lond)
                distda[s] = getpoint.get_distance_grcircle(lata,lona,latd,lond)

                for l in listdictvalues.keys():
                    a = values[indices[0]]
                    b = values[indices[1]]
                    c = values[indices[2]]
                    d = values[indices[3]]
                    print indices,a,b,c,d

        
                    # bilinear:
                    #pab = a * xb[s] / distab[s] + b * xa[s] / distab[s]  
                    #pcd = d * xc[s] / distab[s] + c * xd[s] / distcd[s] 
                    
                    pvalue = (yd[s] * xb[s] / (distab[s] * distda[s]) + 
                              yc[s] * xa[s] / (distab[s] * distbc[s]) + 
                              yb[s] * xd[s] / (distcd[s] * distbc[s]) + 
                              ya[s] * xc[s] / (distcd[s] * distda[s]) )
                    #pvalue = pab * ycd / distda + pcd * yab / distda[s]
                    print 'dus ',a,b,c,d,pvalue
                    print 'en  ',xa[s],xb[s],xc[s],xd[s],xa[s]+xb[s],xc[s]+xd[s]


            if 0: # debug
                print a,b,xb,xa, distab
                print lata,latb,ya,yb,yab
                print latc,latd,yc,yd,ycd
                print pab,pcd,pvalue

            return pvalue

        else: # return average (?)
            return (a+b+c+d)/4.
        
if __name__ == '__main__':

    import pygrib
    import matplotlib.pyplot as plt
    import get_closest_gridpoint  as getpoint

    gribfile = '/home/plas/data/radarver/fc2010071400+017grib'
    gribfile = '/usr/people/plas/python/tools/ex2012062100+012grib'


    #grbs = pygrib.open('test_hir.grb')
    grbs = pygrib.open(gribfile)

    print('===== The grib message considered: =====')
    print grbs[1]
    lats,lons = grbs[1].latlons()
    values  = grbs[1].values
    values2 = grbs[2].values
    values3 = grbs[3].values

    values  = grbs[5].values
    values2 = grbs[6].values
    values3 = grbs[7].values

    grbs.close()

    latt,lont = 52.376,4.879
    latt,lont = 52.372,4.85
    latt,lont = 52.361,4.845
    llstat    = (latt,lont)

    g        = getpoint.grid_searcher(lats,lons,method=2)
    s        = g.get_closest_point(latt,lont)
    a,b,c,d  = g.get_four_surrounding_points(latt,lont)
    pts      = (a,b,c,d)
    indic    = ((pts[0],pts[2]),(pts[0],pts[3]),(pts[1],pts[3]),(pts[1],pts[2]))

    print('\n===== The index and lat,lon of closest gridpoint: =====')
    print latt,lont
    print 'single point: ',s, lats[s],lons[s]
    print 'four points:  ',a,b,c,d

    print('\n===== The indices and values of the surrounding points: =====')
    print indic
    for i in indic:
        print 'lat,lon,value: ',lats[i],lons[i],values[i]

    
    v = interp(values,(lats,lons),llstat,indic)
    print('\n===== The interpolated value: =====')
    print v

    print '_______ using scipy _______'

    if 1:
        x,y,z = [lons[i] for i in indic],[lats[i] for i in indic],[values[i] for i in indic]
        #z = [1,1,2,2]
        xt,yt = lont,latt
        #srbf = Rbf(lats.flatten(),lons.flatten(),values.flatten())
        srbf = Rbf(x,y,z)
        print 'Rbf: ',srbf(xt,yt),' from ',[values[i] for i in indic]
        
        zi = griddata((x, y), z, (xt, yt), method='linear') #'cubic')
        print 'griddata: ',zi,' from ',[values[i] for i in indic]
        
        # define grid
        print 'x',lons[a,c],lons[b,c]
        print 'y',lats[a,c],lats[b,c]
        
        xi = np.linspace(lons[a,c],lons[b,c],10)
        yi = np.linspace(lats[a,c],lats[b,c],10)
        
        xi = np.linspace(4.845,4.88,10)
        yi = np.linspace(52.36,52.38,10)

        #xi = np.linspace(4.842,4.843,10)
        #yi = np.linspace(52.36,52.37,10)

        # grid the data.
        zi = griddata((x,y),z,(xi[None,:], yi[:,None])) 
        print 'gridded, plotting'

        # contour the gridded data, plotting dots at the randomly spaced data points.
        CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
        CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        # plot data points.
        plt.scatter(x,y,marker='o',c='b',s=5)
        plt.scatter(xt,yt,marker='o',c='purple',s=25)
        plt.xlim(min(x),max(x))
        plt.ylim(min(y),max(y))
        plt.title('griddata test') # (%d points)' % npts)
        plt.show()

    # everything:
    if 1:
        print '\n===== using all points ====='
        x,y,z = lats.flatten(),lons.flatten(),values.flatten()

        i1,i2 = 300,300
        j1,j2 = 500,500
        x,y,z = lons[i1:j1,i2:j2].flatten(),lats[i1:j1,i2:j2].flatten(),values[i1:j1,i2:j2].flatten()
        
        xi = np.linspace(lons[i1,i2],lons[i1,j2],100)
        yi = np.linspace(lats[i1,i2],lats[j1,i2],100)
        zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear') #'cubic')
        #zi = np.zeros([10,10])

        CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
        CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar

        xl = [4.879,4.85,4.845]
        yl = [52.376,52.372,52.361]

        plt.scatter(xl,yl,marker='o',c='purple',s=25)
        plt.show()
 

        zi1 = griddata((x, y), z, (xl, yl), method='linear') 
        print 'griddata: ',zi1,' from all points'
        #x,y,z = lats.flatten(),lons.flatten(),values.flatten()

    #sys.exit(0)

    if 0:
        print 10*'=', 'test class: ',10*'='
        latt2,lont2 = 51.9,4.6
        latt3,lont3 = 51.6,4.5
        stlist = ((latt,lont),(latt2,lont2),(latt3,lont3))
        myline = reinterpol(gridll = (lats,lons),latlonlist = stlist)
        myline.compute_dist()

        somevals = {1:values,2:values2,3:values3}
        myline.multi_intp(somevals)
