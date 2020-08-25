#! /usr/bin/env python

import os,sys
import datetime as dt
import numpy as np

#import argparse, 
import pygrib
import netCDF4
import time
from netCDF4 import Dataset, date2num, num2date

from scipy.interpolate import griddata

#from harmoniefields import level_types,field_codes_harmonie
from harmonie_codes import level_types,field_codes_harmonie36,field_codes_harmonie37

# some data, always useful:
stationlist = {
    1:{'name':'Cabauw',
       'sid':  6348,
       'latlon':(51.971, 4.927),
       },
    2:{'name':'Chilbolton',
       'sid': 6348,
       'latlon':(51.1445, -1.437),
       },
    3:{'name':'Palaiseau',
       'sid': 6348,
       'latlon':(48.713, 2.204),
       },
    4:{'name':'Lindenberg',
       'sid': 6348,
       'latlon':(52.1, 14.3),
       },
    5:{'name':'Juelich',
       'sid': 6348,
       'latlon':(50.909237, 6.413622),
       },
    6:{'name':'Schiphol',
       'sid': 6240,
       'latlon':(52.3156,4.7903), 
       },
    7:{'name':'De_Bilt',
       'sid': 6260,
       'latlon':(52.099874, 5.176661), 
       },
    }

# some functions, not methods of the ncprofile class:

def toLatlon(latlons,index):

    return [(latlons[0][i],latlons[1][i]) for i in index]

def do_kdtree(combined_x_y_arrays,points):

    '''
    Use a KDTree method to find the closest gridpoint (from stackoverflow)
    '''

    from scipy.spatial import cKDTree
    
    mytree = cKDTree(combined_x_y_arrays)
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
        for i,j in path:
            if 0 < i < dim0-1 and 0 < j < dim1-1:
                llpath = Path(toLatlon(latlons,path))
                if verb: print 'Corner: ',corner,path,'\n',llpath
            else:
                #break
                print 'Point on edge or beyond, return' 
                return None

        # if the path contains the point of interest, return the indices:
        if llpath.contains_point((xlat,xlon)):
            if verb: print 'Succes! ',corner,path,llpath
            return path

    print 10*'*','No path found!',xlat,xlon
    return None

def intValues(field,latlons,indices,(xlat,xlon),verb=True):
    
    '''
    interpolate using SciPy: griddata(), default at points in xlat,xlon
    Uses latlons and some 2D field
    '''
    
    # recall minimal rectangle:

    # avoid keyError:
    dim0,dim1 = latlons[0].shape
    
    i1,i2 = max(0,min([p[0]-1 for p in indices])), max(0,min([p[1]-1 for p in indices]))
    j1,j2 = min(dim0,max([p[0]+1 for p in indices])), min(dim1,max([p[1]+1 for p in indices]))

    x = latlons[1][i1:j1,i2:j2].flatten()
    y = latlons[0][i1:j1,i2:j2].flatten()
    z = field[i1:j1,i2:j2].flatten()
    
    # use scipy's griddata to interpolate:
    intval = griddata((x, y), z, (xlon, xlat), method='linear')
    
    if verb: print 'from',z,'\n mean is ',z.mean(),' interpol gives ',intval

    return intval


def get_levs_ab_latlon_from_grib(grbfile):

    grbs   = pygrib.open(grbfile)
    grb    = grbs[1]
    nlevs = int(grb.numberOfVerticalCoordinateValues/2 -1)
    pv    = grb.pv
    lats,lons = grb.latlons()
    grbs.close()

    return nlevs,pv,(lats,lons)

def init_index(latlons,stationlist):

    '''
    takes latlons (from grib) and the globally defined "stationlist" 
    to find indices of four surrounding points, to interpolate
    '''
    
    # init indices, ONCE
    for s in stationlist.keys():
        print 'find index for ',stationlist[s]['name'] 
        
        xlat,xlon = stationlist[s]['latlon']
        stationlist[s]['ind'] = find_four_neighbours(latlons,(xlat,xlon),verb=True) #False)
        if stationlist[s]['ind'] == None: 
            del(stationlist[s])
            print stationlist 

    return stationlist # unnecessary? Does it update the global dict()?

def get_varname(parameter,level,leveltype,gribtab,tr = 0,model='harmonie',version='3614'):


    if model.lower() == 'harmonie':
        if int(version) < 3712:
            field_codes = field_codes_harmonie36
        else:
            field_codes = field_codes_harmonie37


    if level_types.has_key(leveltype):
        ltype = level_types[leveltype]
    else: # default
        ltype = '105'   

    # if modellevel: search in index for level 001:
    if ltype == '109': 
        ilevel = 1
    # if pressure level: search for level 100 (?)
    #elif ltype == '105' and int(level) in (100,250,300,500,800,850,900,925,950,1000):
    #    print 'pressure level? ',leveltype # == isobaricInhPa
    #    ilevel = 100
    elif ltype == '100': # and int(level) in (100,250,300,500,800,850,900,925,950,1000):
        print 'Yes? pressure level? ',parameter,leveltype,level # == isobaricInhPa
        ilevel = 001
    elif ltype == '008' and int(level) > 1000:
        print 'TOA radiation'
        ilevel = 000
    else:
        ilevel = level


    # construct index:
    index= "{t:03}-{p:03}-{l:03}-{lt}".format(t=gribtab,p=parameter,l=ilevel,lt=ltype)
    if field_codes.has_key(index):
        print index,field_codes[ index ],parameter,leveltype,level
        try:
            (longname, shortname, units, FAname, standardname) = field_codes[ index ]
        except:
            (longname, shortname, units, FAname) = field_codes[ index ]
            standardname = 'unknown'
    else:
        # for now, just take vars that are described in table
        print 'Not in field_codes:',index
        return None,None,None,None

    if shortname is None or shortname.lower() == 'unknown':
        shortname = index
    elif ltype == '109': 
        shortname += '_ml'
    elif ltype == '100': 
        shortname += '_pl'
    elif ltype == '105' and parameter in (181,184,201):
        if tr == 0:
            shortname += '_inst'
        elif tr == 4:
            shortname += '_acc'
            units = 'kg m**-2'


    if longname is None:
        longname = 'Long Name' #grb.longName
    if units is None:
        units = '?'

    print 'From ',parameter,leveltype,level,' distilled ',longname,shortname,units
    return longname,shortname,units,standardname



class ncprofile:

    # used for converting orography to surface geopotential and v.v.
    standard_g = 9.80665
    pressure_levs = [200,300,400,500,600,700,800,850,900,925,950,1000]

    def __init__(self,startdate,outfile,station,stationid):
        self.startdate = startdate
        self.outfile   = outfile
        self.station   = station
        self.stationid = stationid

        self.initfile()

    def initfile(self,startdate=None,outfile=None):

        if startdate == None: startdate = self.startdate
        if outfile   == None: outfile   = self.outfile

        print 'In initfile:',startdate,outfile

        ## set up a netCDF file
        self.nc = Dataset( outfile, 'w', file_format='NETCDF4' ) 
        self.nc.set_fill_off()

        # create dims
        self.dtime = self.nc.createDimension('time', None)

        # create dimension/variable "station"
        slevel  = self.nc.createDimension('station_id', self.stationid)
        vstation = self.nc.createVariable('station','i4',('station_id',)) 
        vstation[:] = self.stationid
        vstation.long_name = 'Station id'
        vstation.cf_role   = 'timeseries_id'
    
        # create vars
        self.vtime = self.nc.createVariable('time','i4',('time',))
        #self.vtime.units = 'hours since 1950-01-01 00:00:00'
        self.vtime.long_name = 'Valid time'
        self.vtime.standard_name = 'time'
        self.vtime.units = 'seconds since 1970-01-01 00:00:00'

        self.nc.description = 'Vertical profile of model {m} at {loc}'
        self.nc.history = 'Created {dmy}'.format(dmy = dt.datetime.today()) 
        self.nc.institution = 'KNMI, The Netherlands' #, Emiel van der Plas 2014'
        self.nc.title = 'Profile extraction from meso-scale NWP models, KNMI, Emiel van der Plas 2014'
        self.nc.source = 'Model {m}, KNMI'
        self.nc.comment = 'Using a python netCDF4 class (nc_profile_class)'

    def initpv(self,nlevs,ab,latlons):

        self.nlevs   = nlevs
        self.plevs   = self.pressure_levs
        self.latlons = latlons

        dlevel = self.nc.createDimension('level',  self.nlevs)
        hlevel = self.nc.createDimension('hlevel', self.nlevs+1)
        plevel = self.nc.createDimension('plevel', len(self.pressure_levs))

        a    = self.nc.createVariable('a','f4',('hlevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
        b    = self.nc.createVariable('b','f4',('hlevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])

        a[:] = ab[:self.nlevs+1 ]
        b[:] = ab[ self.nlevs+1:]

    def create1D(self,shortname,longname,units,standardname):
        
        vgrib = self.nc.createVariable(shortname,'f4',('station_id','time',),
                               zlib=True,
                               fill_value=netCDF4.default_fillvals['f4'])
        vgrib.longname = longname
        vgrib.units    = units
        vgrib.standardname = standardname

    def create2D(self,shortname,longname,units,standardname):
        
        vgrib = self.nc.createVariable(shortname,'f4',('station_id','time','level'),
                               zlib=True,
                               fill_value=netCDF4.default_fillvals['f4'])

        vgrib.longname = longname
        vgrib.units    = units
        vgrib.standardname = standardname

    def create_pl(self,shortname,longname,units,standardname):
        
        vgrib = self.nc.createVariable(shortname,'f4',('station_id','time','plevel'),
                               zlib=True,
                               fill_value=netCDF4.default_fillvals['f4'])

        vgrib.longname = longname
        vgrib.units    = units
        vgrib.standardname = standardname
        
    def add_value_to_var(self,value,step,level,shortname,longname,units,standardname):

        s = self.stationid #['name']

        #print shortname
        if shortname[-3:] == '_ml': 
            if shortname not in self.nc.variables:
                vgrib = self.create2D(shortname,longname,units,standardname)
            vgrib = self.nc.variables[shortname]
            print 'filling 2D ',s,step,level-1
            vgrib[0,step,level-1] = value

        elif shortname[-3:] == '_pl': 
            if shortname not in self.nc.variables:
                vgrib = self.create_pl(shortname,longname,units,standardname)
            vgrib = self.nc.variables[shortname]
            ind = self.pressure_levs.index(level)
            print 'filling 2D pl',shortname,step,level,ind
            vgrib[0,step,ind] = value
                
        else: 
            if shortname not in self.nc.variables:
                vgrib = self.create1D(shortname,longname,units,standardname)
            vgrib = self.nc.variables[shortname]
            #print 'filling 1D ',step #,level-1
            vgrib[0,step] = value

    def writetime(self,starttime,leadtime):

        validdate = starttime + dt.timedelta(hours = leadtime,seconds=30)
        vtime = self.nc.variables['time']
        writetime = int( date2num(validdate,
                                  vtime.units, 
                                  calendar='standard'))
        vtime[leadtime] = writetime


    def calc_height_using_hydrost_eq(self):
        #  #[ calc. height from pressure p
        phi_surface = self.standard_g * self.orography
        self.calc_geopotential(phi_surface)
        self.convert_geopotential_to_height(phi_surface, self.lat)
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
        
        phi_hl  = np.zeros(self.nlevs+1) # allocate a new array
        ln_p_hl = np.log(self.p_half)      # calculate log of pressure
        phi_hl[self.nlevels] = phi_surface    # init the surface value
        
        # step up in the atmosphere (index 0 is on top!)
        self.calc_virtual_temperature()
        for i in range(self.nlevs-1,0,-1):
            phi_hl[i] = phi_hl[i+1] + \
                        R_dry_air*self.t_v[i]*(ln_p_hl[i+1] - ln_p_hl[i])

        # note: p_hl[0] == 0.
        # so: ln_p_hl[0] will become infinite!
        # to prevent nonsense results, substitute it with a very low
        # but finite number
        ln_p_hl[0] = np.log(1.)
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



    def postproc_radiation(self):

        def deacc(ts):
            # subtract previous value, last value is "off" (?)
            return (ts[1:-1] - ts[:-2])/3600.

        lwd = self.nc.variables['lwdrad']
        swd = self.nc.variables['swdrad']
        lwn = self.nc.variables['lwnrad']
        swn = self.nc.variables['swnrad']

        #self.nc.variables['lwdrad'][:] = [(lwd[i] - lwd[i-1])/3600. for i in range(1,49)]
        #self.nc.variables['swdrad'][:] = [(swd[i] - swd[i-1])/3600. for i in range(1,49)]
        #self.nc.variables['lwnrad'][:] = [(lwn[i] - lwn[i-1])/3600. for i in range(1,49)]
        #self.nc.variables['swnrad'][:] = [(swn[i] - swn[i-1])/3600. for i in range(1,49)]

        self.nc.variables['lwdrad'][:] = deacc(lwd)
        self.nc.variables['swdrad'][:] = deacc(swd)
        self.nc.variables['lwnrad'][:] = deacc(lwn)
        self.nc.variables['swnrad'][:] = deacc(swn)

        swu = self.create1D('swurad', 'Surface SW up radiation', 'W m**-2')
        lwu = self.create1D('lwurad', 'Surface LW up radiation', 'W m**-2')
        netrad = self.create1D('netrad','Surface net total radiation', 'W m**-2')
        self.nc.variables['swurad'][:] = deacc(swd[:]-swn[:]) #[(swd[i] - swn[i] - (swd[i-1] - swn[i-1]))/3600. for i in range(1,49)]
        self.nc.variables['lwurad'][:] = deacc(lwd[:]-lwn[:]) #[(lwd[i] - lwn[i] - (lwd[i-1] - lwn[i-1]))/3600. for i in range(1,49)]
        self.nc.variables['netrad'][:] = deacc(lwn[:]+swn[:]) #[(lwn[i] + swn[i] - (lwn[i-1] + swn[i-1]))/3600. for i in range(1,49)]

    def postproc_uv(self):

        def todegree(u,v,offset):
            # offset -1.5 for u>0,v < 0, else 0.5 
            return -180*(np.arctan2(v,u)/(np.pi) + offset)

        def towindspeed(u,v):
            return np.sqrt(u**2 + v**2)

        def towinddir(u,v):
            condlist = [(v<0) & (u < 0),(v>=0) | (u >= 0)]
            choicelist = [todegree(u,v,+0.5),todegree(u,v,-1.5)]
            return np.select(condlist,choicelist)


        u10   = self.nc.variables['u10']
        v10   = self.nc.variables['v10']
        ugust = self.nc.variables['ugust']
        vgust = self.nc.variables['vgust']
        u_ml  = self.nc.variables['u_ml']
        v_ml  = self.nc.variables['v_ml']

        wind      = self.create1D('wind10',    'Windspeed at 10m',      'm s-2',        'wind_speed')
        winddir   = self.create1D('winddir10', 'Wind direction at 10m', 'degrees north','wind_from_direction')
        gust      = self.create1D('gust',      'Windgust at 10m',       'm s-2',        'wind_speed')
        gustdir   = self.create1D('gustdir',   'Gust direction at 10m', 'degrees north','wind_from_direction')
        windml    = self.create2D('wind_ml',   'Windspeed at ml',       'm s-2',        'wind_speed')
        winddirml = self.create2D('winddir_ml','Wind direction at ml',  'degrees north','wind_from_direction')

        self.nc.variables['wind10'][:]    = towindspeed(u10[:],v10[:])
        #condlist = [(v10[:]<0) & (u10[:] < 0),(v10[:]>=0) | (u10[:] >= 0)]
        #choicelist = [todegree(u10[:],v10[:],+0.5),todegree(u10[:],v10[:],-1.5)]
        self.nc.variables['winddir10'][:] = towinddir(u10[:],v10[:])

        self.nc.variables['gust'][:]    = towindspeed(ugust[:],vgust[:])
        self.nc.variables['gustdir'][:] = towinddir(ugust[:],vgust[:])

        self.nc.variables['wind_ml'][:]    = towindspeed(u_ml[:],v_ml[:])
        self.nc.variables['winddir_ml'][:] = towinddir(u_ml[:],v_ml[:])

    def postproc_pz(self):

        # use a,b and surface pressure to calculate p at half levels, then full levels:
        
        p_surf  = self.nc.variables['Psurf'][0,:]
        a,b     = self.nc.variables['a'][:],self.nc.variables['b'][:]

        # give variables same dimension using numpy.tile:
        self.steps = len(p_surf)
        hlevels = len(a)
        
        ## calculate pressure at half, and full levels
        psml    = np.tile(p_surf[:],(hlevels,1)).T
        aml,bml = np.tile(a,(self.steps,1)),np.tile(b,(self.steps,1))
        p_half  = aml + psml * bml
        p_full  = (p_half[:,:-1] + p_half[:,1:])/2.

        # write as nc variable
        p_ml    = self.create2D('p_ml', 'Pressure at ml', 'hPa','air_pressure')
        self.nc.variables['p_ml'][0,:,:] = p_full

        ## calculate thickness/height of levels: 
        t       = self.nc.variables['T_ml'][0,:,:]
        q       = self.nc.variables['qv_ml'][0,:,:]
        tv      = t + (1 + 0.61 * q)   # virtual temperature
        oro     = self.nc.variables['z_oro'][:]
        #ta      = tv[:,1:] + tv[:,:-1]
        
        # value taken from: McIlveen, Fundamentals of weather and Climate, p.74
        R_dry_air = 287. # [J/(K.kg)]
        
        phi_hl  = np.zeros((self.steps,self.nlevs+1)) # allocate a new array
        ln_p_hl = np.log(p_half)        # calculate log of pressure
        phi_hl[:,self.nlevs] = self.standard_g * oro # init the surface value
        
        # step up in the atmosphere (index 0 is on top!)
        for i in range(self.nlevs-1,0,-1):
            phi_hl[:,i] = phi_hl[:,i+1] + R_dry_air * tv[:,i]*(ln_p_hl[:,i+1] - ln_p_hl[:,i])

        # note: p_hl[0] == 0. so: ln_p_hl[0] will become infinite!
        # to prevent nonsense results, substitute it with a very low
        # but finite number
        ln_p_hl[:,0] = np.log(1.)
        phi_hl[:,0]  = phi_hl[:,1] + R_dry_air * tv[:,0]*(ln_p_hl[:,1] - ln_p_hl[:,0])

        # simple, see profile_new.py for more thorough, latitude dependent calculation
        z_hl = phi_hl/self.standard_g 
        z_mid = (z_hl[:,:-1] + z_hl[:,1:])/2
        z_ml    = self.create2D('z_ml', 'Height at ml', 'm','height')
        self.nc.variables['z_ml'][0,:,:] = z_mid


    def postproc_tq(self):

        # dewpoint and  qv 
        # source http://en.wikipedia.org/wiki/Dew_point
 
        t_0 = 273.15
        
        '''
        wim:
        esat = (0.61078*np.exp(17.2694*(t.sdata-273.16)/(t.sdata-35.86)))
	qsat = (0.622 * esat/(p/1000-esat))
	rh = q.sdata/qsat 

        cisco:
        tmelt = 273.15
        ;do i=0,time-1
        ;  do j=0,nstat-1
        ;    tt=t2m(i,j)-tmelt                                     ; tt in [C]
        ;    es2m(i,j)=6.112*exp(17.67*tt/(243.5+tt))              ; es2m in [hPa]
        ;    e2m(i,j) =es2m(i,j)*rh2m(i,j)/100.                    ; e2m in [hPa]
        ;    td2m(i,j)=td+tmelt                                    ; td2m in [K]
        ;    td=(243.5*log(e2m(i,j)/6.112))/(17.67-log(e2m(i,j)/6.112))
        ;    qv2m(i,j)=.622*es2m(i,j)*1e2/ps(i,j)                  ; qv in [kg/kg]
        ;  end do
        ;end do

        '''

        def td(t,q,p):
            tt   = t[:] - t_0

            if 1: # test
                esat = 0.61078 * np.exp(17.2694 * (tt)/(tt + t_0 - 35.86)) # Bij Cisco 6.122 als factor?
                qsat = 0.622 * esat/((1.e-3)*p[:] - esat)
                rh   = q[:]/qsat 
                #td   = t_0 + (243.5*np.log(esat * rh/6.112))/(17.67- np.log(esat * rh/6.112))

                a,b,c = 6.112, 17.67, 243.5
                ps = a * np.exp( b* tt/ (c + tt))
                pa = rh * ps
                pa[pa<=0] = 1.e-10
                td = c * np.log(pa/a)/(b - np.log(pa/a))

            return rh,td

        td_ml    = self.create2D('Td_ml', 'Dew point temperature', 'K','dew_point_temperature')
        rh_ml    = self.create2D('rh_ml', 'Relative humidity', '%','relative_humidity')

        t       = self.nc.variables['T_ml'][:]
        q       = self.nc.variables['qv_ml'][:]
        p       = self.nc.variables['p_ml'][:] # created in postproc_pz

        rh,td   = td(t,q,p)

        self.nc.variables['td_ml'][:] = td + t_0
        self.nc.variables['rh_ml'][:] = rh


        #t2m     = self.nc.variables['t2m'][:]
        #q2m     = self.nc.variables['q2m'][:]
        #td2m     = self.create1D('plas_td2m',  'Dew point temperature', 'K')
        #self.nc.variables['plas_td2m'][:] = td(t2m,q2m,p[:,-1])


    def add_vars_to_nc(self,grbfile):

        '''
        old version, not in use
        '''

        grbs = pygrib.open(grbfile)
        for grb in grbs:

            par        = grb.indicatorOfParameter
            leveltype  = grb.indicatorOfTypeOfLevel
            levelctype = grb.typeOfLevel
            level      = grb.level
            timerange  = grb.timeRangeIndicator
            gribtab    = grb.table2Version
            
            vals       = grb.values
  
            longname,shortname,units = get_varname(par,
                                                   level,
                                                   levelctype,
                                                   gribtab,
                                                   tr=timerange,
                                                   model='harmonie')

            # actual interpolation
            for s in self.stationlist.keys():
                
                xlat,xlon = self.stationlist[s]['latlon']
                iv = intValues(vals,latlons,self.stationlist[s]['ind'],(xlat,xlon),verb=False)
                
                if shortname[:-3] == '_ml': 
                    if shortname not in self.nc.variables:
                        vgrib = self.create2D()
                    vgrib = self.nc.variables[shortname]
                    vgrib[h,level-1] = iv

                elif shortname[:-3] == '_pl': 
                    if shortname not in self.nc.variables:
                        vgrib = self.create_pl()
                    vgrib = self.nc.variables[shortname]
                    vgrib[h,level-1] = iv
                        
                else: 
                    vgrib = self.create1D()

                vgrib = self.nc.variables[shortname]
