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
#  #[ definition of constant parameters
# see for the list of HARMONIE parameter definitions:
# https://hirlam.org/trac/wiki/HarmonieSystemDocumentation/
#         Forecast/Outputlist
# units are given in the HIRLAM list at:
# https://hirlam.org/trac/wiki/HirlamSystemDocumentation/
#         Forecast/Outputlist
#
par_T     =  11 # temperature          (index   5 in the list) [K]
par_phi_s =   6 # surface geopotential (index   9 in the list) [m^2/s^2]
par_P_s   =   1 # surface pressure     (index  20 in the list) [pa]
par_hum   =  51 # specific humidity    (index 160 in the list) [kg/kg]
# par_RH_s  =  52 # surface relative humidity (index 161/162 in the list)
par_u     =  33 # u wind component     (index 163 in the list) [m/s]
par_v     =  34 # v wind component     (index 164 in the list) [m/s]
par_vv    =  40 # vert. velocity       (index 168 in the list) [not given]
par_TKE   =  200 # turbulent kinetic energy 
par_PCP   =   62 # rain
par_GRP   =  201 # graupel
par_SNW   =   79 # snow
par_CLW   =   76 # Cloud water content
par_CLI   =   58 # Cloud ice 
# HIRLAM uses par. 39 in stead:
# pressure coordinate vertical velocity, with unit [Pa/s]
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
            # EvdP BUG!? llpath = Path([(lats[i],lons[i]) for i in path])
            llpath = Path([(self.lats[i],self.lons[i]) for i in path])
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
        # EvdP BUG!? ind = lli.find_four_neighbours((latt,lont)) #,verb=True)
        ind = self.find_four_neighbours((latt,lont)) #,verb=True)
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
        
class press_alt_provider:

    # conversion constant
    deg2rad = numpy.pi/180.

    # used for converting orography to surface geopotential and v.v.
    standard_g = 9.80665

    def __init__(self, gribfile, lat, lon, orography=0.):
        #  #[ prepare the profile
        # copy inputs
        self.lat = lat
        self.lon = lon
        self.orography = orography # relative to the ellipsoid!
        
        # init to missing/None values
        self.nlevels = None
        self.p_half  = None
        self.p_mid   = None
        self.phi_s   = None
        self.P_s     = None
        self.T       = []
        self.hum     = []
        self.t_v     = None

        # do a dummy read just to retrieve the lats,lons and a,b arrays
        dummyvals,lats,lons = get_grib_values(gribfile,par_T,
                                              leveltype='sfc',level=0,TR=0)
        
        # do this just once for a given profile
        lli = LatLonInterpolator(lats, lons)
        
        grbs = pygrib.open(gribfile)
        for grb in grbs:
            if (grb.levelType == 'sfc'):
                # fill surface values
                if (grb.indicatorOfParameter == par_phi_s):
                    self.phi_s = lli.interpolate(grb.values,lat,lon)
                if (grb.indicatorOfParameter == par_P_s):
                    self.P_s = lli.interpolate(grb.values,lat,lon)
            
            if (grb.levelType == 'ml'):
                if self.p_mid is None:
                    # do this just once for a given grib file
                    self.calc_pressure_levels(grb)
                    
                # TODO1: get these lat,lon values from the radiosonde data
                # TODO2: use 2 grib files for time interpolation
                if (grb.indicatorOfParameter == par_T):
                    # EvdP self.T.append(lli.interpolate(grb.values,lat_ml,lon_ml))
                    self.T.append(lli.interpolate(grb.values,self.lat,self.lon))
                if (grb.indicatorOfParameter == par_hum):
                    # EvdP self.hum.append(lli.interpolate(grb.values,lat_ml,lon_ml))
                    self.hum.append(lli.interpolate(grb.values,self.lat,self.lon))

        grbs.close()

        self.T   = numpy.array(self.T)
        self.hum = numpy.array(self.hum)

        self.calc_height_using_hydrost_eq()
        #  #]

    def calc_pressure_levels(self, grb, P_s=101325.):
        #  #[ get pressure levels for HARMONIE model
        # see http://www.ecmwf.int/research/ifsdocs/DYNAMICS/
        #            Chap2_Discretization4.html#961180
        # for these equations
        ab = grb['pv']
        
        nlevels = len(ab)/2-1
        a = ab[:nlevels+1]
        b = ab[nlevels+1:]
        self.p_half = a + b*P_s
        
        a_mid = (a[:-1] + a[1:])/2
        b_mid = (b[:-1] + b[1:])/2
        self.p_mid = a_mid+b_mid*P_s

        self.nlevels = nlevels
        
        # print 'ab = ',ab
        # print 'p_half = ',p_half
        # print 'p_mid  = ',p_mid
        #  #]

    def calc_grav_acc(self, lat):
        #  #[ calc g using a parametrisation
        #--------------------------------------------------------------
        # code borrowed from the subroutine sgps_calc_z()
        # written by  J.C.W. (John) de Vries, KNMI, 08/04/2003
        #--------------------------------------------------------------
        
        # local variables and parameters
        a1 =  5.2885E-3
        a2 = -5.9E-6
        ge =  9.780356 # in [m/s^2]
        
        lat_radians = lat*self.deg2rad 
        
        sinlat      = numpy.sin(lat_radians)
        sinlat2     = sinlat*sinlat
        
        sin2lat     = numpy.sin(2.*lat_radians)
        sin2lat2    = sin2lat*sin2lat
        
        g_surface = ge*(1. + a1*sinlat2 + a2*sin2lat2)
        
        self.g_surface = g_surface # in [m/s^2]
        #  #]

    def calc_earth_radius(self, lat):
        #  #[ calculate earth shape (sea level) using an ellipsoid shape
        #--------------------------------------------------------------
        # code borrowed from the subroutine sgps_calc_z()
        # written by  J.C.W. (John) de Vries, KNMI, 08/04/2003
        #--------------------------------------------------------------
        R_earth_equatorial = 6378.1e3 # [m] ! re
        R_earth_polar      = 6356.8e3 # [m] ! rp
        R_factor_sq        = (R_earth_equatorial/R_earth_polar)**2
        
        lat_radians = lat*self.deg2rad 
        
        sinlat   = numpy.sin(lat_radians)
        sinlat2  = sinlat*sinlat
        
        coslat   = numpy.cos(lat_radians)
        coslat2  = coslat*coslat
        
        r_no_oro = R_earth_equatorial/numpy.sqrt(R_factor_sq*sinlat2+coslat2)
        
        self.r_no_oro = r_no_oro # in [m]
        #  #]

    def calc_virtual_temperature(self):
        #  #[ calculate virtual temperature
        # t   ! input in [K]
        # h   ! input in [kg/kg]
        # t_v ! result (in [K?])
        
        # see syllabus Fysische Meteorologie I, H.R.A.Wessels, KNMI TR-140, p.20
        # (in my collection: map 1, article nr. 14)
        epsilon = 0.62 
        a0 = (1./epsilon) - 1. # = 0.61290
        
        # convert humidity from [g/kg] to relative humidity in [kg/kg]
        # q = h*1.e-3
        # not needed, HARMONIE provides q directly
        q = self.hum
        
        # calculate virtual temperature to correct for the changing
        # air density at constant pressure for different air compositions
        # (this is a small correction, typically less then 1%)
        self.t_v = (1.+a0*q)*self.T
        #  #]

    def calc_geopotential(self, phi_surface):
        #  #[ calc geopotential phi from pressure p
        # remark: fl indicates the full levels, hl the half levels
        # the nr of full levels is nlevels, the nr of half levels is nlevels+1
        # array: p_hl   ! input in [pa]
        # array: t_fl   ! input in [K]
        # array: hum_fl ! input in [kg/kg]
        # phi_surface   ! input in [m^2/s^2]
        #
        # array phi_hl  ! result in [m^2/s^2]
        
        # value taken from: McIlveen, Fundamentals of weather and Climate, p.74
        R_dry_air = 287. # [J/(K.kg)]
        
        phi_hl  = numpy.zeros(self.nlevels+1) # allocate a new array
        ln_p_hl = numpy.log(self.p_half)      # calculate log of pressure
        phi_hl[self.nlevels] = phi_surface    # init the surface value
        
        # step up in the atmosphere (index 0 is on top!)
        self.calc_virtual_temperature()
        for i in range(self.nlevels-1,0,-1):
            phi_hl[i] = phi_hl[i+1] + \
                        R_dry_air*self.t_v[i]*(ln_p_hl[i+1] - ln_p_hl[i])

        # note: p_hl[0] == 0.
        # so: ln_p_hl[0] will become infinite!
        # to prevent nonsense results, substitute it with a very low
        # but finite number
        ln_p_hl[0] = numpy.log(1.)
        phi_hl[0] = phi_hl[1] + \
                    R_dry_air*self.t_v[0]*(ln_p_hl[1] - ln_p_hl[0])
        
        self.phi_hl = phi_hl
    #  #]
    
    def convert_geopotential_to_height(self, phi_surface, lat):
        #  #[ calc. height z from geopotential phi
        # see my earlier fortran work in:
        #    versiebeheer_mercurial/adm/src/CalcHeightFromP/
        #    CalcHeightUsingHydrostEq.F90
        
        # phi         ! input in [m^2/s^2]
        # phi_surface ! input in [m^2/s^2]
        # lat         ! input in [deg]
        # z           ! result in [m]
        
        # this only holds in case the orography is defined 
        # with respect to the ellipsoid !!!!!
        self.calc_earth_radius(lat)
        r_s = self.r_no_oro + self.orography # [m]
        
        # calculate the local gravitational accelleration at ellipsoid level
        # using a simple parametrization. This is then used as surface value 
        # for g.
        self.calc_grav_acc(lat)
        g_s = self.g_surface  # [m/s^2]
        
        # transform phi to z (see my notebook p. 74)
        self.z_hl = r_s*self.phi_hl/(g_s*r_s-self.phi_hl)    # [m]
        
        self.z_mid = (self.z_hl[:-1] + self.z_hl[1:])/2
        #  #]

    def calc_height_using_hydrost_eq(self):
        #  #[ calc. height from pressure p
        phi_surface = self.standard_g*self.orography
        self.calc_geopotential(phi_surface)
        self.convert_geopotential_to_height(phi_surface, self.lat)
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

if __name__ == '__main__':

    ## get data from tar archives using:
    # /net/bhw379/nobackup_1/users/plas/verif/BULL/extract_fc_from_tar.py

    ## get data from current run from /net/bens03/harmonie_data/GVDB/ 
    ## ==>NB files with HARM_N25 are on limited 300x300 domain
    ## at highest resolution of 2.5 x 2.5 km and
    ## contain all fields on high resolution, including DZDT
    ## ==>HARM_N55* files contain a larger area at lower resolution
    ## of 5.5 x 5.5 km
    ## ==>HARM_L25* files contain a larger area at lower resolution
    ## of 5.5 x 5.5 km, but only a selection of a few fields.
    
    harmonie_data_dir = '/net/bhw379/nobackup/users/plas'
    #harmonie_data_dir = '/nobackup/users/kloedej/Harmonie_Data_Examples'
    gribfile = os.path.join(harmonie_data_dir,'fc2013061806+012grib_37h12')
    #gribfile = os.path.join(harmonie_data_dir,'HARM_N25_201310070000_00100_GB')
    
    ml_pars           = [ par_T,   par_hum,   par_u,   par_v,   par_vv]
    ml_par_names      = ['par_T', 'par_hum', 'par_u', 'par_v', 'par_vv']

    # do a dummy read just to retrieve the lats,lons and a,b arrays
    dummyvals,lats,lons = get_grib_values(gribfile,par_T,
                                          leveltype='sfc',level=0,TR=0)
    
    # do this just once for a given profile
    lli = LatLonInterpolator(lats, lons)

    nlevels = 60
    levels = range(1,nlevels+1)

    ml_profile_values = {}
    for ml_par_name in ml_par_names:
        ml_profile_values[ml_par_name] = numpy.zeros(nlevels)

    grbs = pygrib.open(gribfile)
    for grb in grbs:
        if (grb.levelType == 'ml'):
            # TODO1: get these lat,lon values from the radiosonde data
            # TODO2: use 2 grib files for time interpolation
            lat_ml,lon_ml = 52.3,4.2
            for ml_par, ml_par_name in zip(ml_pars,ml_par_names):
                if (grb.indicatorOfParameter == ml_par):
                    #print 'filling level {} for par {}'.\
                    #      format(grb.level, ml_par_name)
                    ml_profile_values[ml_par_name][grb.level-1] = \
                               lli.interpolate(grb.values,lat_ml,lon_ml)

    t_fl   = ml_profile_values['par_T']
    hum_fl = ml_profile_values['par_hum']

    lat_ml,lon_ml = 52.3,4.2
    pa_prov = press_alt_provider(gribfile,lat_ml,lon_ml)

    for i in range(pa_prov.nlevels+1):
        if i<pa_prov.nlevels:
            print "i=",i," p_hl = ",pa_prov.p_half[i],\
                  'phi_hl[i]=',pa_prov.phi_hl[i],\
                  " t_fl = ",pa_prov.T[i]," hum_fl = ",pa_prov.hum[i],\
                  " z_hl[i] = ",pa_prov.z_hl[i]
    
    import matplotlib, matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    new_fontsize = 6
    #            ax    title             xlabel              par
    plotdefs = [(ax1, 'profile of T',   'T [K]',            'par_T'),
                (ax2, 'profile of hum', 'humidity [kg/kg]', 'par_hum'),
                (ax3, 'profile of u',   'u [m/s]',          'par_u'),
                (ax4, 'profile of vv',  'vertical velocity [m/s]', 'par_vv')]
    for ax, title, xlabel, par in plotdefs:
        ax.plot(ml_profile_values[par],pa_prov.z_mid*1.e-3,
                marker='o', markersize=2)
        ax.set_title(title,fontsize=new_fontsize)
        ax.set_xlabel(xlabel,fontsize=new_fontsize)
        ax.set_ylabel('Altitude [km]',fontsize=new_fontsize)
        for l in ax.get_xticklabels()+ax.get_yticklabels():
            l.set_fontsize(fontsize=new_fontsize)
        
    plotfile = 'testprofile.png'
    fig.savefig(plotfile, dpi = 150)
    os.system('xv '+plotfile)
