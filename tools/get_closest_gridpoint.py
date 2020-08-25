#!/usr/bin/env python

#  #[ documentation
"""
A class to easily find the index to the closest grid-point
for a given lat-lon pair, even for a non-rectangular grid.

Assumptions:
-the data is organised in some regular grid, even if they do
 not form straight but curved lines in a lat-lon coordinate system.
-the data may define a global or a local field
-no assumption is made on which array index corresponds to lat or lon
 (since this is not useful for projected or rotated grids anyway)
-input lat,lon arrays and lat,lon points are assumed to be in degrees

Todo's:
-test and implement exceptions for:
 ==>poles (lat=-90 or 90)
 ==>date-boundary (lon=-180 or 180)
 ==>requesting points outside the defined domain
"""
#  #]
#  #[ imported modules
import os, sys, numpy
import pygrib
import matplotlib, matplotlib.pyplot
import unittest
#  #]
#  #[ some constants
dummy_method       = 0
brute_force_method = 1
walker_method      = 2
subdomain_method  = 3

deg2rad = numpy.pi/180.
rad2deg = 180./numpy.pi

R_earth_equatorial = 6378.1
# R_earth_polar      = 6356.8_r8_ ! [km]
# source: http://wwwflag.wr.usgs.gov/USGSFlag/Space/wall/earth.html
#  #]

# a custom exception
class GridSearchException(Exception):
    """ an exception for custom use by the grid_searcher class """
    pass

# this lat,lon distance calculation is copied from the example
# in our genscat and aeolus software
def latlon2xyz_radius1(lat,lon):
    #  #[
    # convert spherical to cartesian coordinates, i.e.
    # calculate x,y,z coordinates for a given lat,lon pair on a sphere
    # of radius 1
    #
    # equations taken from: http://mathworld.wolfram.com/GreatCircle.html

    lon_rad = lon*deg2rad
    lat_rad = lat*deg2rad
    coslat  = numpy.cos(lat_rad)
    x       = numpy.cos(lon_rad)*coslat
    y       = numpy.sin(lon_rad)*coslat
    z       = numpy.sin(lat_rad)
    return (x,y,z)
    #  #]
def get_angle_distance_grcircle(lat1,lon1,lat2,lon2):
    #  #[
    # determine the angle (in degrees) on the earths surface
    # between 2 sets of lat-lon coordinates along a great-circle 
    # (assuming that the shape of the earth is a perfect globe)
    #
    # equations taken from: http://mathworld.wolfram.com/GreatCircle.html
    #
    # use the property that the angle alpha between 2 vectors in
    # cartesian coordinates is given by the inner product of the 2 vectors:
    #              _    _
    # cos(alpha) = v1 . v2
    #
    # note that the error for using a sperical earth in stead of the proper
    # spheroid must always be between the result using R_earth_equatorial
    # and the result using R_earth_polar. Therefore the maximum error
    # is estimated at: 
    # 100.*(R_earth_equatorial - R_earth_polar)/R_earth_polar = 0.34 percent
    # which is an acceptable accuracy for our use

    x1,y1,z1 = latlon2xyz_radius1(lat1,lon1)
    x2,y2,z2 = latlon2xyz_radius1(lat2,lon2)
    inner_prod = (x1*x2+y1*y2+z1*z2)

    # extra check, abs(inner_prod) should never be larger than 1 !
    # however, due to rounding errors it might happen anyway
    if (abs(inner_prod) > 1.): inner_prod = 1.

    # extra check, abs(inner_prod) should never be smaller than -1 !
    # however, due to rounding errors it might happen anyway
    if (abs(inner_prod) < -1.): inner_prod = -1.

    # note that an inner product may be negative
    # so apply abs() to ensure the distance is always positive or zero
    angle_dist = abs(numpy.arccos(inner_prod)*rad2deg)
    return angle_dist
    #  #]
def get_distance_grcircle(lat1,lon1,lat2,lon2):
    #  #[
    # this function converts the angle-distance between 
    # 2 latlon pairs to a km distance
    # since angle_distance will always yield a non-negative
    # number, distance will always be non-negative as well
    angle_distance = get_angle_distance_grcircle(lat1,lon1,lat2,lon2) 
    distance = deg2rad*R_earth_equatorial*angle_distance
    return distance
    #  #]

class grid_searcher:
    def __init__(self,lats,lons,method=brute_force_method,
                 allow_outside_domain=False,globalfield=False,
                 verbose=False):
        #  #[
        # some sanity checks on the inputs
        if lats.shape!=lons.shape:
            txt = 'ERROR in __init__: '+\
                  'lats and lons arrays should have identical shape!'
            raise GridSearchException(txt)
        if method not in [dummy_method, brute_force_method,
                          walker_method, subdomain_method]:
            txt = 'ERROR in __init__: '+\
                  'unknown method!'
            raise GridSearchException(txt)

        # store the inputs
        self.lats   = lats
        self.lons   = lons
        self.method = method
        self.allow_outside_domain = allow_outside_domain
        self.verbose = verbose
        
        # some derived values
        self.npoints_axis1 = lats.shape[0]
        self.npoints_axis2 = lats.shape[1]
        self.num_distances_calculated = 0

        self.globalfield = globalfield
        self.wraplons = False        
        if (globalfield):
            # assume we have a global field!
            self.wraplons = True
            # in this case we expect a non-rotated regular lon grid
            # ranging either from 0 ... 360-dlon
            # or ranging from -180 ... 180-dlon
            # for both cases, add an extra column to self.lats and
            # self.lons, and correct j coordinate afterwards
            self.lats = numpy.zeros((lats.shape[0],lats.shape[1]+1))
            self.lons = numpy.zeros((lons.shape[0],lons.shape[1]+1))
            self.lats[:,:lats.shape[1]] = lats
            self.lons[:,:lons.shape[1]] = lons
            self.lats[:,lats.shape[1]] = self.lats[:,0]
            self.lons[:,lats.shape[1]] = self.lons[:,0]+360.
            self.npoints_axis2 = lats.shape[1]+1
            self.lon_range = [numpy.min(self.lons),
                              numpy.max(self.lons)]
        #  #]
    def get_four_surrounding_points(self,lat,lon):
        #  #[
        i,j = self.get_closest_point(lat,lon)
        # try the 4 boxes around the found location to see
        # in which one the requested lat,lon is located

        #debug = True
        debug = False

        #  #[ lower left
        ibox = [i-1,i]
        jbox = [j-1,j]
        if debug:
            print 'debug: ======================================='
            print 'debug: trying: lower left; i=',i,'j=',j
            print 'debug: trying: ibox=',ibox,'jbox=',jbox
        if self.wraplons:
            # take care of lon wrapping
            if j==0 or j==self.npoints_axis2-1:
                jbox = [self.npoints_axis2-2,0]
                if debug:
                    print 'debug: take care of lon wrapping'
                    print 'debug: modified: jbox=',jbox
                    
        if ibox[0]>=0 and jbox[0]>=0:
            if self.inside_box(lat,lon,ibox,jbox,debug=debug):
                if debug: print 'debug: using lower left box'
                return ibox[0],ibox[1],jbox[0],jbox[1]
            else:
                if debug: print 'NOT INSIDE'
        #  #]

        #  #[ lower right
        ibox = [i-1,i]
        jbox = [j,j+1]
        if debug:
            print 'debug: ======================================='
            print 'debug: trying: lower right; i=',i,'j=',j
            print 'debug: trying: ibox=',ibox,'jbox=',jbox
        if self.wraplons:
            # take care of extra row of lon values
            if j == self.npoints_axis2-2:
                jbox = [self.npoints_axis2-2,0]
                if debug:
                    print 'debug: take care of extra row of lon values'
                    print 'debug: modified: jbox=',jbox
            
        if ibox[0]>=0 and jbox[1]<self.npoints_axis2-1:
            if self.inside_box(lat,lon,ibox,jbox,debug=debug):
                if debug: print 'debug: using lower right box'
                return ibox[0],ibox[1],jbox[0],jbox[1]
            else:
                if debug: print 'NOT INSIDE'
        #  #]

        #  #[ upper left
        ibox = [i,i+1]
        jbox = [j-1,j]
        if debug:
            print 'debug: ======================================='
            print 'debug: trying: upper left; i=',i,'j=',j
            print 'debug: trying: ibox=',ibox,'jbox=',jbox
            
        if self.wraplons:
            # take care of lon wrapping
            if j==0 or j==self.npoints_axis2-1:
                
                jbox = [self.npoints_axis2-2,0]
                if debug:
                    print 'debug: take care of lon wrapping'
                    print 'debug: modified: jbox=',jbox

        if ibox[1]<self.npoints_axis1-1 and jbox[0]>=0:
            if self.inside_box(lat,lon,ibox,jbox,debug=debug):
                if debug: print 'debug: using upper left box'
                return ibox[0],ibox[1],jbox[0],jbox[1]
            else:
                if debug: print 'NOT INSIDE'
        #  #]

        #  #[ upper right
        ibox = [i,i+1]
        jbox = [j,j+1]
        if debug:
            print 'debug: ======================================='
            print 'debug: trying: upper right; i=',i,'j=',j
            print 'debug: trying: ibox=',ibox,'jbox=',jbox
            
        if self.wraplons:
            # take care of extra row of lon values
            if j == self.npoints_axis2-2:
                jbox = [self.npoints_axis2-2,0]
                if debug:
                    print 'debug: take care of extra row of lon values'
                    print 'debug: modified: jbox=',jbox

        if ibox[1]<self.npoints_axis1-1 and jbox[1]<self.npoints_axis2-1:
            if self.inside_box(lat,lon,ibox,jbox,debug=debug):
                if debug: print 'debug: using upper right box'
                return ibox[0],ibox[1],jbox[0],jbox[1]
            else:
                if debug: print 'NOT INSIDE'
        #  #]

        txt = 'ERROR: this point should not be reached'
        if debug:
            print 'debug: i,j = ',i,j
            print 'debug: ibox = ',ibox
            print 'debug: jbox = ',jbox
            res = self.inside_box(lat,lon,ibox,jbox,debug=True)
            print 'debug: res = ',res       
            print txt
            os._exit(1)
            
        raise GridSearchException(txt)
        #  #]
    def get_closest_point(self,lat,lon):
        #  #[
        if self.method == brute_force_method:
            i,j = self.get_closest_point_brute_force(lat,lon)
        elif self.method == dummy_method:
            i,j = self.get_closest_point_dummy(lat,lon)
        elif self.method == walker_method:
            i,j = self.get_closest_point_walk(lat,lon)
        elif self.method == subdomain_method:
            i,j = self.get_closest_point_subdomains(lat,lon)

        #if self.wraplons:
        #    # take care of extra row of lon values
        #    print 'wrapping around ! i,j=',i,j
        #    print 'self.npoints_axis2 = ',self.npoints_axis2
        #    
        #    if j == self.npoints_axis2-1: j=0
            
        # check if the requested point is inside the defined latlon domain
        if not self.allow_outside_domain:
            #print 'debug, allow_outside_domain:'
            #print 'i,j = ',i,j
            ibox = [i-1,i+1]
            jbox = [j-1,j+1]
            if i==0:
                ibox=[i,i+1]
            if i==self.npoints_axis1-1:
                ibox=[i-1,i]
            if j==0:
                jbox=[j,j+1]
            if j==self.npoints_axis2-1:
                jbox=[j-1,j]

            if not self.inside_box(lat,lon,ibox,jbox):
                print 'debug: ibox,jbox = ',ibox,jbox
                txt = 'ERROR in get_closest_point: '+\
                      'requested latlon seems outside domain!'
                raise GridSearchException(txt)

        # finally compensate for the extra lon-row if needed
        if self.wraplons:
            if j==self.npoints_axis2-1:
                j=0
        return i,j
        #  #]
    def get_closest_point_dummy(self,lat,lon):
        return (0,0)
    def get_distance_i_j(self,lat,lon,i,j):
        #  #[
        if i<0 or i > self.npoints_axis1-1:
            return None
        if j<0 or j > self.npoints_axis2-1:
            return None
        lat2 = self.lats[i,j]
        lon2 = self.lons[i,j]
        distance = get_distance_grcircle(lat,lon,lat2,lon2)
        self.num_distances_calculated += 1
        return distance
        #  #]
    def get_closest_point_brute_force(self,lat,lon,ibox=None,jbox=None):
        #  #[ just walk along all points to find the closest
        closest_distance = None
        closest_ii = None
        closest_jj = None
        if ibox:
            istart, iend = ibox
        else:
            istart, iend = 0, self.npoints_axis1-1
        if jbox:
            jstart, jend = jbox            
        else:
            jstart, jend = 0, self.npoints_axis2-1
        for ii in range(istart,iend+1):
            for jj in range(jstart,jend+1):
                distance = self.get_distance_i_j(lat,lon,ii,jj)
                if closest_distance is None:
                    closest_distance = distance
                    closest_ii = ii
                    closest_jj = jj
                else:
                    if distance<closest_distance:
                        #print 'closer point found; distance: ',distance,\
                        #      'ii,jj=',ii,jj
                        closest_distance = distance
                        closest_ii = ii
                        closest_jj = jj

        return (closest_ii,closest_jj)
        #  #]
    def get_closest_point_walk(self,lat,lon):
        #  #[ walk along x or y axis depending on which is closest

        # start in the center of the grid
        i = int(self.npoints_axis1/2.)
        j = int(self.npoints_axis2/2.)
        
        progressing = True

        while progressing:
            # check left, right, up and down (if possible)
            distance       = self.get_distance_i_j(lat,lon,i,j)
            distance_left  = self.get_distance_i_j(lat,lon,i-1,j  )
            distance_right = self.get_distance_i_j(lat,lon,i+1,j  )
            distance_up    = self.get_distance_i_j(lat,lon,i  ,j+1)
            distance_down  = self.get_distance_i_j(lat,lon,i  ,j-1)
            # print 'distance       = ',distance
            # print 'distance_left  = ',distance_left
            # print 'distance_right = ',distance_right
            # print 'distance_up    = ',distance_up
            # print 'distance_down  = ',distance_down
            if   distance_left  and (distance_left <distance): i=i-1
            elif distance_right and (distance_right<distance): i=i+1
            elif distance_up    and (distance_up   <distance): j=j+1
            elif distance_down  and (distance_down <distance): j=j-1
            else: progressing = False
            #
            # print 'i,j,dist=',i,j,distance
        return (i,j)
        #  #]
    def get_closest_point_subdomains(self,lat,lon,subdomain=None,prefix=''):
        #  #[ split in 4 boxes and choose recursively

        if subdomain:
            domain = subdomain
        else:
            domain = (0,0,self.npoints_axis1-1,self.npoints_axis2-1)

        if self.verbose:
            print prefix+'trying domain: ',domain

        if ( (domain[2]-domain[0]<=3) or 
             (domain[3]-domain[1]<=3)   ):
            if self.verbose:
                print 'i or j size too small, stopping iteration'
            return self.get_closest_point_brute_force(lat,lon,\
                                    ibox=(domain[0],domain[2]),
                                    jbox=(domain[1],domain[3]))
        
        icenter = domain[0]+int((domain[2]-domain[0])/2.)
        jcenter = domain[1]+int((domain[3]-domain[1])/2.)

        domain_upper_left  = (domain[0],domain[1],  icenter,  jcenter)
        domain_lower_left  = (domain[0],  jcenter,  icenter,domain[3])
        domain_upper_right = (  icenter,domain[1],domain[2],  jcenter)
        domain_lower_right = (  icenter,  jcenter,domain[2],domain[3])

        for (domain_to_try,name) in [(domain_upper_left,'domain_upper_left'),
                                     (domain_lower_left,'domain_lower_left'),
                                     (domain_upper_right,'domain_upper_right'),
                                     (domain_lower_right,'domain_lower_right')]:
            if self.verbose:
                print prefix+'sub: ',name,' trying domain: ',domain_to_try
            is_inside = self.inside_box(lat,lon,
                                        (domain_to_try[0],domain_to_try[2]),
                                        (domain_to_try[1],domain_to_try[3]),
                                        prefix=prefix+'==>')
            if is_inside:
                # recursively call myself untill convergence
                if self.verbose:
                    print prefix+'INSIDE!'
                return self.get_closest_point_subdomains(\
                            lat,lon,subdomain=domain_to_try,
                            prefix=prefix+'  ')
            else:
                if self.verbose:
                    print prefix+'NOT inside domain: ',domain_to_try
                
        # requested point seems to be outside valid domain!
        #
        txt = 'ERROR in get_closest_point_subdomains: '+\
              'requested latlon seems outside domain!'
        raise GridSearchException(txt)
    
        #  #]
    def inside_box(self,lat,lon,ibox,jbox,prefix='',debug=False):
        #  #[
        # a simplified method, may be inaccurate for very
        # distorted boxes !!!
        #print 'debug: ibox, jbox = ',ibox,jbox
        minlat = min( self.lats[ibox[0],jbox[0]],
                      self.lats[ibox[0],jbox[1]],
                      self.lats[ibox[1],jbox[0]],
                      self.lats[ibox[1],jbox[1]] )
        maxlat = max( self.lats[ibox[0],jbox[0]],
                      self.lats[ibox[0],jbox[1]],
                      self.lats[ibox[1],jbox[0]],
                      self.lats[ibox[1],jbox[1]] )
        minlon = min( self.lons[ibox[0],jbox[0]],
                      self.lons[ibox[0],jbox[1]],
                      self.lons[ibox[1],jbox[0]],
                      self.lons[ibox[1],jbox[1]] )
        maxlon = max( self.lons[ibox[0],jbox[0]],
                      self.lons[ibox[0],jbox[1]],
                      self.lons[ibox[1],jbox[0]],
                      self.lons[ibox[1],jbox[1]] )

        wrapped_lon = lon
        if self.wraplons:
            if debug:
                print 'debug: lat,lon = ',lat,lon
                print 'debug: minlat,maxlat = ',minlat,maxlat
                print 'debug: minlon,maxlon = ',minlon,maxlon
                print 'debug: self.lon_range = ',self.lon_range

            # take care of lon being outside -180..180 or 0..360 domain
            if lon < self.lon_range[0]:
                wrapped_lon = lon+360.
            if lon > self.lon_range[1]:
                wrapped_lon = lon-360.
            if debug:
                print 'debug: wrapped_lon=',wrapped_lon

            # take care of swapped lon indices
            if (maxlon - minlon) > 270.:
                tmp    = maxlon
                maxlon = minlon + 360.
                minlon = tmp
                if debug:
                    print 'debug: swapped min/max lon'
                    print 'debug: minlon,maxlon = ',minlon,maxlon
            
        if self.verbose:
            print prefix+'lat = ',lat,'minlat,maxlat = ',minlat,maxlat
            print prefix+'lon = ',lon,'minlon,maxlon = ',minlon,maxlon
            print prefix+'wrapped_lon = ',wrapped_lon

        outside_box = (lat<minlat or lat>maxlat or
                       wrapped_lon<minlon or wrapped_lon>maxlon)
        return (not outside_box)
        #  #]
        
if __name__ == '__main__':
    #  #[ default settings+test definitions
    do_plots       = False
    do_manualtests = False
    do_unittests   = False

    print_num_distances = False

    # method_to_use = brute_force_method # takes 87 s.
    # method_to_use = walker_method # takes 0.53 s.
    method_to_use = subdomain_method # takes 0.011s!!!

    testfiles     = ['lambert.grb','reglatlon.grb','rotll.grb',
                     'global_reglatlon.grb']

    # pairs of (lat,lon) in deg.
    testlocations_lambert   = [ [(50.,2.),( 52.,-4.)], # inside domain
                                [(58.,-2.),(48.,0.)] ] # outside domain
    testlocations_reglatlon = [ [(50.,2.), (52.,6.)],  # inside domain
                                [(58.,-2.),(48.,3.)] ] # outside domain
    testlocations_rotll     = [ [(50.,2.), (70.,-56.)],  # inside domain
                                [(20.,0.),(30.,75.)] ] # outside domain
    testlocations_global_reglatlon = \
          [ [(50.,2.), (30.,75.),(10.,359.9),             
             (-10.,-35.),(-25.,-179.75),(89.9, 31.),
             (-71.7,0.05),
             (10.1,359.9),(9.9,359.9),(10.1,0.1),(9.9,0.1),
             (10.1,-0.1),(9.9,-0.1),
             (-55.2,-0.37),
             (54.4, 0.32),
             ],  # inside domain
            [(91.,0.),(-95.,75.)] ] # outside domain
    
    testlocations = [testlocations_lambert,
                     testlocations_reglatlon,
                     testlocations_rotll,
                     testlocations_global_reglatlon]
    
    #testlocations = [ 
    #                  [(50.,2.),(30.,75.),(10.,359.9),(-10.,-35.),
    #                   (-25.,-179.75), (89.9, 31.), ]]
    
    # expected index and closest points # i,j,i1,j1,i2,j2
    expected_test_results_lambert   = [ (( 60, 171),(59,170,60,171)),
                                        (( 151, 4), (150,4,151,5)) ]
    expected_test_results_reglatlon = [ (( 43, 56),  (43,55,44,56)),
                                        (( 130, 167),(130,166,131,167)) ]
    expected_test_results_rotll     = [ (( 147, 341),(146,341,147,342)),
                                        (( 332, 217),(331,216,332,217)) ]
    expected_test_results_global_reglatlon = \
                        [ (( 280, 4),(279,3,280,4)),
                          (( 240, 150),(239,149,240,150)),
                          (( 200, 0),  (199,719,200,0)),
                          (( 160, 650),(159,649,160,650)),
                          (( 130, 361),(129,360,130,361)),
                          (( 360, 61), (359,61,360,62)),
                          ((37,0),(36,0,37,1)),
                          ((200,0),(200,719,201,0)),
                          ((200,0),(199,719,200,0)),
                          ((200,0),(200,0,201,1)),
                          ((200,0),(199,0,200,1)),
                          ((200,0),(200,719,201,0)),
                          ((200,0),(199,719,200,0)),
                          ((70,719),(69,719,70,0)),
                          ((289,1),(288,0,289,1)),
                        ]
    expected_test_results = [expected_test_results_lambert,
                             expected_test_results_reglatlon,
                             expected_test_results_rotll,
                             expected_test_results_global_reglatlon]
    #  #]
    #  #[ commandline handling

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if arg=='plot':
                do_plots = True
            elif arg=='test':
                do_unittests = True
            elif arg=='manual':
                do_manualtests = True
            else:
                print 'Undefined option: ',arg
                print 'Usage: '
                print sys.argv[0]+' [plot|manual|test]'
                sys.exit(1)

    #  #]
    for (fnr,testfile) in enumerate(testfiles):
        #  #[ prepare test
        print '+'*75
        print 'loading grib file: ',testfile
        print '+'*75
        # open a GRIB file, create a grib message iterator:
        grbs = pygrib.open(testfile)
        grb =  grbs.next()
        # get the latitudes and longitudes of the grid:
        lats, lons = grb.latlons()
        print 'lats.shape = ',lats.shape
        print 'lons.shape = ',lons.shape
        #  #]
        if do_plots:
            #  #[ do some test plots
            print '+'*75
            # prepare plotting a test plot
            fig = matplotlib.pyplot.figure()
            ax = fig.add_subplot(1,1,1) # rows,columns,count
            ax.plot(lons[::10,::10],lats[::10,::10],
                    marker='.',markersize=1.5,linestyle='')
            plotfile = 'testplot'+testfile+'.png'
            sf = fig.savefig(plotfile)
            print 'created plot: '+plotfile
            os.system('xv '+plotfile+'&')
            #  #]
        if do_manualtests:
            #  #[ do some manual tests
            print '+'*75
            print 'doing manual tests for testfile: ',testfile
            print 'testlocations = ',testlocations[fnr]
            globalfield = testfile == 'global_reglatlon.grb'
            g = grid_searcher(lats,lons,method=method_to_use,
                              globalfield=globalfield)
            for (lat,lon) in testlocations[fnr][0]:
                print '='*50
                print 'trying to find lat,lon: ',lat,lon,' for grid: ',testfile
                print '='*50
                try:
                    i,j = g.get_closest_point(lat,lon)
                    print 'i,j = ',i,j
                    print 'lat,lon = ',lat,lon
                    print 'closest gridpoint is: ',\
                          lats[i,j],lons[i,j]
                    print 'nr. of distances that have been calculated: ',\
                          g.num_distances_calculated
                except GridSearchException as e:
                    print 'GridSearchException(1) detected ...'
                    #if fnr==3:
                    #    print 'running again to display proper traceback'
                    #    i,j = g.get_closest_point(lat,lon)
                    #    #raise e
                try:
                    i1,i2,j1,j2 = g.get_four_surrounding_points(lat,lon)
                    print 'i1,i2,j1,j2 = ',i1,i2,j1,j2
                    print 'box boundaries: '
                    print 'a:',lats[i1,j1],lons[i1,j1]
                    print 'b:',lats[i1,j2],lons[i1,j2]
                    print 'c:',lats[i2,j1],lons[i2,j1]
                    print 'd:',lats[i2,j2],lons[i2,j2]
                except GridSearchException as e:
                    print 'GridSearchException(2) detected ...'
                    #if fnr==3:
                    #    print 'running again to display proper traceback'
                    #    i1,i2,j1,j2 = g.get_four_surrounding_points(lat,lon)
                    raise e
            #  #]
        if do_unittests:
            #  #[ do some unit tests
            print '+'*75
            print 'setting up unittests for testfile: ',testfile

            #class CheckGridSearcher(unittest.TestCase):
            #    def test_inside_domain(self):
            #        g = grid_searcher(lats,lons,
            #                          method=brute_force_method)
            #        lat,lon = testlocations[fnr][0]
            #        i,j = g.get_closest_point(lat,lon)
            #        n = g.num_distances_calculated
            #        self.assertEqual((i,j,n),
            #                         expected_test_results[fnr][0])

            # the following code does the same, but automatically
            # constructs 3 different classes with 3 different names
            # in the global namespace, and adds individual testfunctions
            # to them.
            # see: http://www.pythonexamples.org/2011/01/12/\
            #      how-to-dynamically-create-a-class-at-runtime-in-python/

            # note: defining optional parameters fnr=fnr and friends
            # is needed here
            # to ensure these values are saved! Otherwise each of the
            # test functions will use the values defined in the last run
            # of the for loop and will run 4 times the same test!
            def test_inside_domain(self,fnr=fnr,lats=lats,lons=lons,
                                   testfile=testfile):
                #print 'testing file nr. ',fnr
                globalfield = testfile == 'global_reglatlon.grb'
                g = grid_searcher(lats,lons,method=method_to_use,
                                  globalfield=globalfield)
                ncases = len(testlocations[fnr][0])
                for i in range(ncases):
                    (lat,lon) = testlocations[fnr][0][i]
                    (exp_i,exp_j), dummy = expected_test_results[fnr][i]

                    i,j = g.get_closest_point(lat,lon)
                    if print_num_distances:
                        print 'num_distances_calculated = ',\
                              g.num_distances_calculated
                    self.assertEqual((i,j),(exp_i,exp_j))

            def test_outside_domain(self,fnr=fnr,lats=lats,lons=lons,
                                    testfile=testfile):
                globalfield = testfile == 'global_reglatlon.grb'
                g = grid_searcher(lats,lons,method=method_to_use,
                                  globalfield=globalfield)
                ncases = len(testlocations[fnr][1])
                for i in range(ncases):
                    (lat,lon) = testlocations[fnr][1][i]
                    self.assertRaises(GridSearchException,
                                      g.get_closest_point,lat,lon)

            def test_get_4_surrounding_points(self,fnr=fnr,lats=lats,lons=lons,
                                             testfile=testfile):
                globalfield = testfile == 'global_reglatlon.grb'
                g = grid_searcher(lats,lons,method=method_to_use,
                                  globalfield=globalfield)

                ncases = len(testlocations[fnr][0])
                for i in range(ncases):
                    (lat,lon) = testlocations[fnr][0][i]
                    dummy, exp_box = expected_test_results[fnr][i]
                    
                    i1,i2,j1,j2 = g.get_four_surrounding_points(lat,lon)
                    self.assertEqual((i1,j1,i2,j2),exp_box)
                
#    except GridSearchException as e:
#                print 'GridSearchException detected ...'

            base,ext = os.path.splitext(testfile)
            CheckGridSearcher = \
                    type('CheckGridSearcher_'+str(base),
                         (unittest.TestCase,),
                         {'test_inside_domain':test_inside_domain,
                          'test_outside_domain':test_outside_domain,
                          'test_get_4_surrounding_points':
                          test_get_4_surrounding_points}
                         )

            globals()['CheckGridSearcher_'+str(base)] = CheckGridSearcher

            # needed to prevent the last one appearing as an additional
            # double test case
            del(CheckGridSearcher)
            #  #]
        #  #[ close down test
        grbs.close()
        #  #]

    if do_unittests:
        # needed because my own argument handling confuses the
        # unittest system
        sys.argv=[sys.argv[0],]
        
        # run all tests
        print 'running unittests:'
        unittest.main()
