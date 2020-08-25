#! /usr/bin/env python

#  #[ documentation
# First version written by Emiel van der Plas, KNMI, 2-Oct-2013
# at my request.
# Adapted by Jos de Kloe for requesting model properties collocated
# with an actual radiosonde trajectory that varies in location and time.
#  #]
#  #[ imported modules
import os,sys
import numpy
import scipy.interpolate, scipy.spatial
# Path is used to determine if point is inside path: (uses pnpoly, presumably)
from matplotlib.path import Path
import pygrib
#  #]

class LatLonInterpolator:
    def __init__(self, lats, lons, distance_upper_bound=0.2, verbose=False):
        #  #[ init the class
        # limit search distance to assure no artefacts occur in case
        # of rotated grids that may look curved in a rectangular latlon
        # representation

        # store the inputs
        self.distance_upper_bound = distance_upper_bound
        self.lats = lats
        self.lons = lons
        self.verbose = verbose
        # hardcoded for now
        self.modelname='harmonie'
        #self.modelname='ECMWF'
        
        # store shape of lat and lon arrays
        self.dim0, self.dim1 = self.lats.shape
        self.ngridpoints = self.dim0 * self.dim1
        # print 'dim0,dim1,dim0*dim1 = ',self.dim0, self.dim1, self.ngridpoints

        # flatten the lat and lon arrays
        self.lats_1d = self.lats.ravel()
        self.lons_1d = self.lons.ravel()

        # combine multiple 1D arrays into a single 2D array
        llarray = numpy.dstack([self.lats_1d, self.lons_1d])[0]

        # prepare a KDTree to allow finding the closest gridpoint
        self.mytree = scipy.spatial.cKDTree(llarray)
        #  #]
    def convert_index_1d_to_2d(self,i):
        #  #[ convert index in 1D arrays into index pairs in 2D arrays
        i0,j0 = i % self.dim1, (i-i % self.dim1)/self.dim1
        if self.modelname=='harmonie':
            return j0,i0
        else: # case for Herlam/D11, ECMWF
            return i0,j0
        #  #]
    def find_four_neighbours(self,(xlat,xlon),verb=True):
        #  #[ find neighbours
        '''
        Find the indices of the four surrounding gridpoints to do the
        interpolation using the values at these corners
        Uses the lats,lons and a lat,lon of a point of interest
        '''

        # query the KDTree to find the closest gridpoint
        dist, res2 = self.mytree.query((xlat,xlon),
                          distance_upper_bound=self.distance_upper_bound)
        # print 'res2 = ',res2, 'dist = ',dist

        if res2 == self.ngridpoints:
            print 'Sorry, no closest point or neighbours found within the'
            print 'constraint: distance_upper_bound = {}'.\
                  format(self.distance_upper_bound)
            return None
        
        # convert index in 1D arrays into index pairs in 2D arrays
        i0,j0 = self.convert_index_1d_to_2d(res2)
    
        if not ( (self.lats_1d[res2] == self.lats[i0,j0]) and
                 (self.lons_1d[res2] == self.lons[i0,j0])     ):
            print 'Sorry, index definition is wrong, '+\
                  'maybe you used a wrong modeltype?'
            print 'KDTree: i0,j0,res2 = ',i0,j0,res2
            print 'lat/lon 1d:',self.lats_1d[res2],self.lons_1d[res2]
            print 'lat/lon 2d:',self.lats[i0,j0],self.lons[i0,j0]
            print 'These lat/lon pairs should be identical!'
            return None
        
        # determine if point is inside path: (uses pnpoly, presumably)
        # from matplotlib.path import Path
        
        # define 4 paths on all sides of nearest neighbour:
        sq = {}
        sq['ll']   = [(i0,  j0  ),(i0+1,j0  ),(i0+1,j0+1),(i0,  j0+1)]
        sq['lr']   = [(i0-1,j0  ),(i0,  j0  ),(i0,  j0+1),(i0-1,j0+1)]
        sq['ur']   = [(i0-1,j0-1),(i0,  j0-1),(i0,  j0  ),(i0-1,j0  )]
        sq['ul']   = [(i0,  j0-1),(i0+1,j0-1),(i0+1,j0  ),(i0,  j0  )]
        
        # loop over 4 adjacent grid cells:
        for corner in sq:
            path = sq[corner]
            llpath = Path([(lats[i],lons[i]) for i in path])
            if self.verbose: 
                print 'For point ',xlat,xlon, 'from nn ',i0,j0
                print 'Corner: ',corner,path,'\n',llpath

            # if the path contains the point of interest, return the indices:
            if llpath.contains_point((xlat,xlon)):
                if self.verbose: 
                    print 'Succes! ',corner,path,llpath
                    # for ll in list(llpath):
                    #     print 'So backsubstitute in lats,lons ',\
                    #           ll,lats[ll],lons[ll]
                return path

        print 10*'*','No path found!',xlat,xlon
        return None
        #  #]
    def interpolate(self,values,latt,lont):
        #  #[ do the actual interpolation
        ind = lli.find_four_neighbours((latt,lont)) #,verb=True)
        if ind is not None:
            # print 'found neighbour index pairs:', ind
            pass
        else:
            print 'no neighbour index pairs found'
            sys.exit(1)

        # instead of taking the whole field to scipy.griddata, we take the
        # four surrounding points
        # if the field supplied to griddata is too large, a memory leak is
        # there to kill your job
        
        # copy the four corner values of the surrounding box
        # to a small list before calling interpolate
        vals  = [values[ij] for ij in ind]
        # print vals
        
        zlats = [self.lats[ij] for ij in ind]
        zlons = [self.lons[ij] for ij in ind]
        
        dmo = scipy.interpolate.griddata((zlons,zlats),vals,
                                         (lont,latt),method='linear')
        #print str((dmo, grb.dataDate,grb.dataTime,grb.level))#+'\r',
        #sys.stdout.flush()
        return dmo
        #  #]
        
def get_grib_values(gribfile,par=61,leveltype='sfc',level=457,TR=0):
    #  #[ retrieve grib fields
    grbs = pygrib.open(gribfile)
    for grb in grbs:
        #print grb
        #print grb.indicatorOfParameter,grb.levelType, \
        #      grb.level,grb.timeRangeIndicator 
        if ( (grb.indicatorOfParameter == par) and
             (grb.levelType == leveltype)      and
             (grb.level == level)              and
             (grb.timeRangeIndicator == TR)        ):
            # print 'found field: ',grb
            vals = grb.values
            lats,lons = grb.latlons()
            break # first one should be the one

    grbs.close()
    return vals,lats,lons
    #  #]

def get_pressure_levels(grb, P_s=101325.):
    #  #[ get pressure levels for HARMONIE model
    ab = grb['pv']
    nlevels = len(ab)/2-1
    a = ab[:nlevels+1]
    b = ab[nlevels+1:]
    #for aa,bb in zip(a,b):
    #    print '{:21.15f} {:21.15f}'.format(aa,bb)
    a_mid = (a[:-1] + a[1:])/2
    b_mid = (b[:-1] + b[1:])/2
    p_mid = a_mid+b_mid*P_s
    return p_mid
    #  #]
    
if __name__ == '__main__':

    ## get data from tar archives using:
    # /net/bhw379/nobackup_1/users/plas/verif/BULL/extract_fc_from_tar.py

    ## get data from current run from /net/bens03/harmonie_data/GVDB/ 
    ## NB files with HARM_N25 are on limited 300x300 domain
    ## but contain all fields on high resolution, including DZDT

    harmonie_data_dir = '/net/bhw379/nobackup/users/plas'
    #harmonie_data_dir = '/nobackup/users/kloedej/Harmonie_Data_Examples'
    gribfile = os.path.join(harmonie_data_dir,'fc2013061806+012grib_37h12')

    par = 11 # T
    #par = 40 # VV in _md files

    # do a dummy read just to retrieve the lats,lons and a,b arrays
    dummyvals,lats,lons = get_grib_values(gribfile,par=181,
                                          leveltype='sfc',level=0,TR=0)
    
    # do this just once for a given profile
    lli = LatLonInterpolator(lats, lons, distance_upper_bound=0.2)

    grbs = pygrib.open(gribfile)
    nlevels = 60
    levels = range(1,nlevels+1)
    profile_values = numpy.zeros(nlevels)
    p_array = None
    for grb in grbs:
        if ( (grb.indicatorOfParameter == par) and
             (grb.levelType == 'ml')           and
             (grb.level in levels)                 ):
            if p_array is None:
                # do this just once for a given grib file
                p_array = get_pressure_levels(grb)
            # TODO1: get these lat,lon values from the radiosonde data
            # TODO2: use 2 grib files for time interpolation
            latt,lont = 52.3,4.2
            profile_values[grb.level-1] = lli.interpolate(grb.values,latt,lont)

    import matplotlib, matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax1 = fig.add_subplot(111)
    #ax2 = ax1.twinx()
    ax1.set_title('profile of T')
    ax1.set_xlabel('T [K]')
    ax1.set_ylabel('Pressure [hPa]')
    #ax1.semilogy(profile_values,p_array*1.e-2)
    ax1.plot(profile_values,p_array*1.e-2)
    ax1.set_yscale('log')
    ax1.set_ylim(1100.,1.)
    #ax1.set_ylim(max(p_array),min(p_array))
    plotfile = 'testprofile.png'
    fig.savefig(plotfile) # , dpi = 150)
    os.system('eog '+plotfile)
