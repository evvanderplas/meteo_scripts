#! /usr/bin/env python

import os,sys,glob
import datetime as dt
import numpy as np
from os.path import join as pjoin 

#import argparse, 
import pygrib
import netCDF4
import time
from netCDF4 import Dataset #, date2num, num2date

from scipy.interpolate import griddata

#from harmoniefields import level_types,field_codes_harmonie
from harmonie_codes_kpt import level_types,field_codes_harmonie36,field_codes_harmonie37

# debugging:
from memory_inspector import report_mem_usage

# some data, always useful:
stationlist = {
    1:{'name':'Cabauw',
       'sid':  6348,
       'latlon':(51.971, 4.927),
       'z':0.81,
       },
    2:{'name':'Chilbolton',
       'sid': 6348,
       'latlon':(51.1445, -1.437),
       'z':45.,
       },
    3:{'name':'Palaiseau',
       'sid': 7157,
       'latlon':(48.713, 2.204),
       'z':112.,
       },
    4:{'name':'Lindenberg',
       'sid': 10393,
       'latlon':(52.1, 14.3),
       'z':103.81,
       },
    5:{'name':'Juelich',
       'sid': 10503,
       'latlon':(50.909237, 6.413622),
       'z':82.,
       },
    6:{'name':'Schiphol',
       'sid': 6240,
       'latlon':(52.3156,4.7903), 
       'z':-5.,
       },
    7:{'name':'De_Bilt',
       'sid': 6260,
       'latlon':(52.099874, 5.176661), 
       'z':5.,
       },
    8:{'name':'P11B',
       'sid': 6203,
       'latlon':(52.3517, 3.3347), 
       'z':0.,
       },
    }

# some functions, not methods of the ncprofile class:

def toLatlon(latlons,index):

    return [(latlons[0][i],latlons[1][i]) for i in index]

def toposix(dtobj,starttime=dt.datetime(1970,1,1,0,0,0)):

    if type(dtobj) != dt.datetime:
        print 'something wrong with datetime:',dtobj
        sys.exit(1)

    timedl = dtobj - starttime

    return int(timedl.total_seconds())


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

def get_varname(parameter,level,leveltype,gribtab,tr = 0,model='harmonie',version='3614',verb=False):


    if model.lower() == 'harmonie':
        if int(version) < 3712 or gribtab == 1:
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
    elif ltype == '100': # and int(level) in (100,250,300,500,800,850,900,925,950,1000):
        #print 'Yes? pressure level? ',parameter,leveltype,level # == isobaricInhPa
        ilevel = 001
    elif ltype == '008' and int(level) > 1000:
        print 'TOA radiation'
        ilevel = 000
    else:
        ilevel = level

    # construct index:
    index= "{t:03}-{p:03}-{l:03}-{lt}".format(t=gribtab,p=parameter,l=ilevel,lt=ltype)
    ## Harmonie code table does not (yet) specify time range indicator: 
    ## for now: add special here, later: redraw harmonie_codes table!
    if parameter in (181,184,201) and ltype == 105 and tr == 0:
        index= "{t:03}-{p:03}-{l:03}-{lt}-{tr}".format(t=gribtab,p=parameter,l=ilevel,lt=ltype,tr=tr)
        
    ##############################################
    # Parse time length indicator for Harmonie37
    ##############################################
    
    sprefix,lprefix,lsuffix = "","",""
    if int(version) > 3711: # new gribtable, harmonie_codes_37
        #tri = grb.timeRangeIndicator
        if tr == 0:
            lprefix = 'Instantaneous '
        elif tr == 2:
            lsuffix = ' over interval'
        elif tr == 4:
            sprefix = 'a'
            lprefix = 'Accumulated '
        else:
            print "Cannot parse timeRangeIndicator ", tr
    
    if field_codes.has_key(index):
        if verb: print index,field_codes[ index ],parameter,leveltype,level
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
    if longname is None:
        longname = 'Long Name' #grb.longName
    if units is None:
        units = '?'

    shortname = sprefix+shortname
    longname  = lprefix+longname+lsuffix

    if verb: print 'From ',parameter,leveltype,level,' distilled ',longname,shortname,units
    return longname,shortname,units,standardname

class ncconv:

    # used for converting orography to surface geopotential and v.v.
    standard_g = 9.80665
    pressure_levs = [200,300,400,500,600,700,800,850,900,925,950,1000]

    def __init__(self,startdate,indir = None,infiles = [],info=None,outdir = './',stationlist = stationlist):

        # gather some information
        self.startdate = startdate
        self.ymdh = self.startdate.strftime('%Y%m%d%H')
        self.indir     = indir
        self.outdir    = outdir
        self.stationlist = stationlist
        self.infiles   = infiles
        self.info      = info


        try:
            os.makedirs(self.outdir)
        except OSError:
            print 'Outdir exists: ',self.outdir, os.path.isdir(self.outdir)



    def init_info(self,force = False):

        '''
        Figuring out how files are prsented to module
        either using only the info (dict) with or without a preformatted directory, 
        or presenting files.

        Not yet fool-proof!
        '''

        #outfilepattern = 'profile_{model}_{ymdh}_{loc}.nc'.format(model=self.info['model'],ymdh= self.ymdh,loc = '*')
        outfilepattern = 'allvar_{loc}_{ymdh}_{model}.nc'.format(model=self.info['model'],ymdh= self.ymdh,loc = '*')
        if len(glob.glob(pjoin(self.outdir,outfilepattern))) > 0 and not force: 
            print 'Profiles: Already done: back to you Jim!' #; sys.exit(0)
            return 'Done'

        if self.infiles == []:
            if self.indir == None and self.info != None:
                y,m,d,h   = self.startdate.strftime('%Y'),self.startdate.strftime('%m'),self.startdate.strftime('%d'),self.startdate.strftime('%H')
                self.indir = self.info['dir'].format(ymdh = self.ymdh, dtgh=self.ymdh, y=y, m=m, d=d, h=h)
            if self.info != None:
                fformat,ltmax = self.info['fformat'],self.info['ltmax']
                lastfile = pjoin(self.indir,fformat.format(ymdh = self.ymdh,dtgh = self.ymdh,lt = str(ltmax).zfill(3)))
        else:
            if info == None:
                self.info = {'model':'default',
                             'version':3811,
                             'name':'Harmonie',
                             'ltmax':48,
                         }
            self.infiles = infiles
            lastfile = self.infiles[-1]

        # 
        if not os.path.exists(lastfile): 
            print 'No file',lastfile # /net/bxsfsn_stc/nfs_stc/lam_oper/prodharm/scratch/hm_home/38h11/archive/{y}/{m}/{d}/{h}/HA38_N25_201409040000_04800_GB
            return 'Done' #'File definition awry' # sys.exit(1)
        
        self.nlevs,self.ab,self.latlons = get_levs_ab_latlon_from_grib(lastfile)
        self.stationlist                = init_index(self.latlons,self.stationlist)

        # a data structure
        self.times  = np.zeros(self.info['ltmax'] + 1, dtype = 'int')

        # write time anyway
        for lt in range(0,self.info['ltmax']+1,1):
            self.times[lt] = toposix(self.startdate + dt.timedelta(hours = lt))

        self.names  = {}
        self.values = {}
        for s in self.stationlist.keys():
            self.values[s] = {}

        # every gribfile? Every leadtime modulo step:
        self.step = 1

        return 'todo'

    def loop_over_grib(self):
        
        for lt in range(0,self.info['ltmax']+1,self.step):
        #for lt in range(0,10): #NB EvdP debug!

            #print 'In nc_class: Looping over grib files'
            #report_mem_usage() # report

            gfile = pjoin(self.indir,self.info['fformat'].format(ymdh = self.ymdh,dtgh = self.ymdh,lt = str(lt).zfill(3)))
            print 'Parsing ',gfile, os.path.exists(gfile),dt.datetime.today()
            if os.path.exists(gfile):
                self.parse_gribfile(gfile,lt)
            else:
                pass

    def parse_gribfile(self,gfile,lt):

        model,version = self.info['model'],self.info['version']
        
        grbs = pygrib.open(gfile)
        for grb in grbs:
        
            par        = grb.indicatorOfParameter
            leveltype  = grb.indicatorOfTypeOfLevel
            levelctype = grb.typeOfLevel
            level      = grb.level
            gribtab    = grb.table2Version
            TR         = grb.timeRangeIndicator

            vals       = grb.values

            #print par,level,levelctype,gribtab,model,version
            #if levelctype == 'ml': print 10*'*', 'ML !'

            longname,shortname,units,stdname = get_varname(par,
                                                           level,
                                                           levelctype,
                                                           gribtab,
                                                           model=model,
                                                           version=version,
                                                           tr = TR
                                                       )

            if shortname == None:
                print 'Again: not in field codes',par,level,levelctype,gribtab
                continue

            for s in self.stationlist.keys():
                xlat,xlon = self.stationlist[s]['latlon']
                ind       = self.stationlist[s]['ind']
                iv        = intValues(vals,self.latlons,ind,(xlat,xlon),verb=False)

                parindex = shortname
                if not self.values[s].has_key(parindex):
                    # get names
                    self.names[parindex] = [longname,shortname,units,stdname]

                    # create array per station, init per type:
                    if levelctype == 'hybrid':
                        self.values[s][parindex] = np.zeros((self.info['ltmax']+1,self.nlevs))
                    elif levelctype == 'isobaricInhPa':
                        self.values[s][parindex] = np.zeros((self.info['ltmax']+1,len(self.pressure_levs)))
                    else:
                        self.values[s][parindex] = np.zeros(self.info['ltmax']+1)

                if levelctype == 'hybrid':
                    self.values[s][parindex][lt,level-1] = iv
                elif levelctype == 'isobaricInhPa':
                    self.values[s][parindex][lt,self.pressure_levs.index(level)] = iv
                else:
                    self.values[s][parindex][lt] = iv

        grbs.close()

    def initprofiles(self,startdate=None,outfile=None):

        #if startdate == None: startdate = self.startdate
        #if outfile   == None: outfile   = self.outfile

        print 'In initfile:',self.startdate #,outfile

        self.files = {}
        ## set up netCDF files
        for s in self.stationlist.keys():
            #outfile = pjoin(self.outdir,'profile_{model}_{ymdh}_{loc}.nc'.format(model=self.info['model'],ymdh= self.ymdh,loc = self.stationlist[s]['name']))
            outfile = pjoin(self.outdir,'allvar_{loc}_{ymdh}_{model}.nc'.\
                                format(model=self.info['model'], 
                                       ymdh=self.ymdh, 
                                       loc=self.stationlist[s]['name'])
                            )

            if os.path.exists(outfile):
                print 'removing old ',outfile
                os.remove(outfile)
            print 'Creating',outfile
            self.files[s] = Dataset( outfile, 'w', file_format='NETCDF4' ) 
            self.files[s].set_fill_off()

            # create dims
            dtime    = self.files[s].createDimension('time', None)
            #dstation = self.files[s].createDimension('station', None)

            # create dimension/variable "station"
            slevel   = self.files[s].createDimension('station', 1) #self.stationid)

            if 0:
                vstation = self.files[s].createVariable('station','str',('station',)) 
                #vstation = self.files[s].createVariable('station','i4',('station',)) 
                vstation[0] = str(self.stationlist[s]['sid'])
                vstation.long_name = 'Station id'
                vstation.cf_role   = 'timeseries_id'
            
            vlat = self.files[s].createVariable('lat','f4',('station',)) 
            vlat[:] = self.stationlist[s]['latlon'][0]
            vlat.long_name = 'Station latitude'
            vlat.units     = 'degrees_north'

            vlon = self.files[s].createVariable('lon','f4',('station',)) 
            vlon[:] = self.stationlist[s]['latlon'][1]
            vlon.long_name = 'Station longitude'
            vlon.units     = 'degrees_east'

            vheight = self.files[s].createVariable('height','f4',('station',)) 
            vheight[:] = self.stationlist[s]['z']
            vheight.long_name = 'Station altitude'
            vheight.units     = 'm'

            rtime = self.files[s].createVariable('forecast_reference_time','i4',())
            rtime.standard_name = 'time'
            rtime.units = 'seconds since 1970-01-01 00:00:00'
            rtime.calendar = 'standard'
            rtime.long_name     = 'forecast_reference_time'
            rtime.standard_name = 'forecast_reference_time'
            rtime[:] = toposix(self.startdate)

            # create vars
            vtime = self.files[s].createVariable('time','i4',('time',))
            vtime.long_name = 'Valid time'
            vtime.standard_name = 'time'
            vtime.units = 'seconds since 1970-01-01 00:00:00'
            vtime[:] = self.times

            self.files[s].description = 'Vertical profile of model {m} at {loc}'.format(m=self.info['model'], loc=self.stationlist[s]['name'])

            gribfiles = self.info['fformat'].format(ymdh = self.ymdh,dtgh = self.ymdh,lt = '*')
            self.files[s].history = 'Created from forecast (GRIB {gribfile}) at {dmy}'.\
                                    format(dmy = dt.datetime.today(), gribfile = gribfiles) 
            self.files[s].institution = 'KNMI, The Netherlands' #, Emiel van der Plas 2014'
            self.files[s].title = 'Harmonie profile extraction, KNMI, Emiel van der Plas 2014'
            self.files[s].source = 'Model {m}, KNMI'.format(m=self.info['model'])
            self.files[s].comment = 'Using a custom python netCDF4 class (nc_class)'
            self.files[s].featureType = 'timeSeries'

            # create dimension levels: #def initpv(self):

            #dlevel = self.files[s].createDimension('level',  self.nlevs)
            #hlevel = self.files[s].createDimension('hlevel', self.nlevs+1)
            # voor KPT:
            dlevel = self.files[s].createDimension('nlev',  self.nlevs)
            hlevel = self.files[s].createDimension('nlevp1', self.nlevs+1)

            plevel = self.files[s].createDimension('plevel', len(self.pressure_levs))
            plevs  = self.files[s].createVariable('plevel','i4',('plevel',))
            plevs.long_name = 'Pressure levels'
            plevs.standard_name = 'pressure'
            plevs.units = 'hPa'
            plevs[:] = self.pressure_levs

            #a    = self.files[s].createVariable('a','f4',('hlevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
            #b    = self.files[s].createVariable('b','f4',('hlevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
            a    = self.files[s].createVariable('a','f4',('nlevp1',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
            b    = self.files[s].createVariable('b','f4',('nlevp1',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
            
            a[:] = self.ab[:self.nlevs+1 ]
            b[:] = self.ab[ self.nlevs+1:]


    def writefile(self):

        # write the values gathered in self.values dict to file
        # access to netCDF file can be slow, so we do it all at once!

        for s in self.stationlist.keys():
            for shortname in self.values[s].keys():
                vals = self.values[s][shortname]
                #print shortname,vals.shape,(self.info['ltmax']+1,),(self.info['ltmax']+1,self.nlevs),(self.info['ltmax']+1,len(self.pressure_levs))
                if vals.shape == (self.info['ltmax']+1,):
                    self.create1D(s,self.files[s],shortname,vals)
                elif vals.shape == (self.info['ltmax']+1,self.nlevs):
                    self.create2D(s,self.files[s],shortname,vals)
                elif vals.shape == (self.info['ltmax']+1,len(self.pressure_levs)):
                    self.create_pl(s,self.files[s],shortname,vals)
                else:
                    print 'Huh?!@%',shortname,vals.shape
                    sys.exit(1)

    def create1D(self,station_index,ncfile,shortname,values):

        netcdf_file = self.files[station_index]
        station_id  = self.stationlist[station_index]['sid']

        #vgrib = self.nc.createVariable(shortname,'f4',('station_id','time',),
        vgrib = netcdf_file.createVariable(shortname,'f4',('station','time',),
                                           zlib=True,
                                           fill_value=netCDF4.default_fillvals['f4'])

        (longname,shortname,units,stdname) = self.names[shortname]

        print 'Check shapes: nc var, values: ', vgrib.shape, values.shape, values.reshape([1]+list(values.shape)).shape
        if vgrib.shape != values.shape:
            values = values.reshape(vgrib.shape)

        vgrib.long_name = longname
        vgrib.units     = units
        vgrib.standardname = stdname
        #vgrib[0,:] = values
        vgrib[:] = values

    def create2D(self, station_index, ncfile, shortname, values): #shortname,longname,units,standardname):
        
        netcdf_file = self.files[station_index]
        station_id  = self.stationlist[station_index]['sid']

        #vgrib = self.nc.createVariable(shortname,'f4',('station_id','time','level'),
        #vgrib = ncfile.createVariable(shortname,'f4',('time','level'), 
        vgrib = netcdf_file.createVariable(shortname,'f4',('station','time','nlev'),
                                           zlib=True,
                                           fill_value=netCDF4.default_fillvals['f4'])

        (longname,shortname,units,stdname) = self.names[shortname]

        if vgrib.shape != values.shape:
            values = values.reshape(vgrib.shape)

        vgrib.long_name = longname
        vgrib.units    = units
        vgrib.standardname = stdname
        #vgrib[0,:] = values
        vgrib[:] = values

    def create_pl(self,station_index,ncfile,shortname,values): #shortname,longname,units,standardname):
        
        #vgrib = self.nc.createVariable(shortname,'f4',('station_id','time','plevel'),
        netcdf_file = self.files[station_index]
        station_id  = str(self.stationlist[station_index]['sid'])

        vgrib = netcdf_file.createVariable(shortname, 'f4', ('station','time','plevel'),
                                      zlib=True, fill_value=netCDF4.default_fillvals['f4'])

        print 'Check shapes: nc var, values: ', vgrib.shape, values.shape

        (longname,shortname,units,stdname) = self.names[shortname]

        if vgrib.shape != values.shape:
            values = values.reshape(vgrib.shape)

        vgrib.long_name = longname
        vgrib.units    = units
        vgrib.standardname = stdname
        #vgrib[str(station_id),:] = values
        vgrib[:] = values #.reshape([1]+list(values.shape))
        
    def writetime(self,starttime,leadtime):

        for f in self.ncfiles:
            vtime    = f.variables['time']
            vtime[:] = self.times


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
        
        delta   = 1.e-10
        phi_hl  = np.zeros(self.nlevs+1) # allocate a new array
        ln_p_hl = np.log(self.p_half + delta) # calculate log of pressure
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
            res = np.zeros(len(ts),dtype=float)
            res[1:-1] = (ts[1:-1] - ts[:-2])/3600.
            #res[0]    = ts[0]
            res[0]    = res[1]
            return res
            
        self.names['swsu']     = ('Surface SW up radiation', 'swsu', 'W m**-2','surface_net_upward_shortwave_flux')
        self.names['lwsu']     = ('Surface LW up radiation', 'lwsu', 'W m**-2','surface_net_upward_longwave_flux')
        self.names['radsd']    = ('Surface net total radiation', 'radsd', 'W m**-2','surface_net_downward_radiative_flux')
        self.names['lhflux']   = ('Latent heat flux through evaporation', 'lhflux', 'W m**-2','latent_heat_flux')

        '''
        error message:

            self.create1D(self.files[s],'SWnet',deacc(swn))
            Post-processing
            /usr/people/plas/python/tools/nc_class.py:801: RuntimeWarning: divide by zero encountered in log
              ln_p_hl = np.log(p_half)        # calculate log of pressure
            Traceback (most recent call last):
              File "/usr/people/plas/python/tools/use_ncprof_cli.py", line 33, in <module>
                nobj.postproc_radiation()
              File "/usr/people/plas/python/tools/nc_class.py", line 686, in postproc_radiation
                self.names['SWnet']    = self.names['aswsn']
            KeyError: 'aswsn'
            Done: Profile loop over leadtimes

        '''

        if 'aswsn' in self.names.keys():
            # timeRangeIndicator was set
            accdict = {'swsn':'aswsn',
                       'lwsn':'alwsn',
                       'swsd':'aswsd',
                       'lwsd':'alwsd',
                       'lhevap':'alhevap',
                       'latf':'alatf',
                       'senf':'asenf',
                       'ustress':'austress',
                       'vstress':'avstress',
                       'waterevap':'awaterevap',
                   }
        else:
            # timeRangeIndicator was NOT set
            accdict = {'swsn':'swsn',
                       'lwsn':'lwsn',
                       'swsd':'swsd',
                       'lwsd':'lwsd',
                       'lhevap':'lhevap',
                       'latf':'latf',
                       'senf':'senf',
                       'ustress':'ustress',
                       'vstress':'vstress',
                       'waterevap':'waterevap',
                   }


        # cut off "Accumulated", (longname,shortname,units,stdname) = self.names[shortname]
        #for var in ['SWnet','LWnet','aswsd','alwsd']: #,'aswsu','alwsu']:
        for var in accdict.keys():
            print var, var in self.names.keys(),accdict[var],accdict[var] in self.names.values()
            #if accdict[var] in self.names.values():
            try:
                longname = self.names[accdict[var]][0] #= self.names[accdict[var]][0][12:]
                print 'From ',longname,' to ',' '.join(longname.split()[1:])
                self.names[accdict[var]][0] = ' '.join(longname.split()[1:])
            except KeyError:
                print 'Name not in parsed grib: ',var

        
        self.names['SWnet']    = self.names[accdict['swsn']]
        self.names['LWnet']    = self.names[accdict['lwsn']]

        for s in self.stationlist.keys():

            for var in ['latf','senf','ustress','vstress','lhevap','waterevap']:
                try:
                    print 'Deaccumulating ',accdict[var], accdict[var] in self.names.keys()
                    print 'Shape orig',self.files[s].variables[accdict[var]].shape
                    new = deacc( self.files[s].variables[accdict[var]][0,:])
                    print 'Shape new',new.shape,new.reshape(self.files[s].variables[accdict[var]].shape).shape
                    self.files[s].variables[accdict[var]][:] =  deacc( self.files[s].variables[accdict[var]][0,:]).reshape(self.files[s].variables[accdict[var]].shape)
                    self.files[s].variables[accdict[var]].long_name = self.names[accdict[var]][0]
                except KeyError:
                    print 'No deaccumulation: not in parsed grib: ',var

            lwd = self.files[s].variables[accdict['lwsd']][0,:] #['lwdrad']
            swd = self.files[s].variables[accdict['swsd']][0,:] #['swdrad']
            lwn = self.files[s].variables[accdict['lwsn']][0,:] #['lwnrad']
            swn = self.files[s].variables[accdict['swsn']][0,:] #['swnrad']

            swnet,lwnet   =  deacc(swn),deacc(lwn)
            swdown,lwdown =  deacc(swd),deacc(lwd)
 
            self.files[s].variables[accdict['swsd']][:] = swdown.reshape(self.files[s].variables[accdict['swsd']].shape)
            self.files[s].variables[accdict['lwsd']][:] = lwdown.reshape(self.files[s].variables[accdict['lwsd']].shape)
            self.create1D(s,self.files[s],'SWnet',swnet)
            self.create1D(s,self.files[s],'LWnet',lwnet)

            self.create1D(s,self.files[s],'swsu' ,swdown - swnet) #deacc(swd[:]-swn[:]))
            self.create1D(s,self.files[s],'lwsu' ,lwdown - lwnet) #deacc(lwd[:]-lwn[:]))
            self.create1D(s,self.files[s],'radsd',lwnet + swnet ) #deacc(lwn[:]+swn[:]))

            if 1: #for var in ['latf','senf','ustress','vstress','lhevap','waterevap']:
                try:
                    self.create1D(s,self.files[s],'lhflux',deacc(self.files[s].variables[accdict['lhevap']][0,:]))
                    #self.create1D(self.files[s],accdict[var],
                except KeyError:
                    print 'No latent heat flux through evaporation'

    def postproc_uv(self):

        def todegree(u,v,offset):
            # offset -1.5 for u>0,v < 0, else 0.5 
            return -180*(np.arctan2(v,u)/(np.pi) + offset)

        def towindspeed(u,v):
            return np.sqrt(u**2 + v**2)

        def towinddir(u,v):
            if 0:
                condlist = [(v<0) & (u < 0),(v>=0) | (u >= 0)]
                choicelist = [todegree(u,v,+0.5),todegree(u,v,-1.5)]
                return np.select(condlist,choicelist)
            else:
                return np.mod(-90.- (180./np.pi) * np.arctan2(v,u),360.) 

        if 0: # names I like but not compatible to testbed?
            self.names['wind10']     = ('Windspeed at 10m',      'wind10',    'm s-2',         'wind_speed') #(longname,shortname,units,stdname)
            self.names['winddir10']  = ('Wind direction at 10m', 'winddir10', 'degrees north', 'wind_from_direction') #(longname,shortname,units,stdname)
            self.names['wind_ml']    = ('Windspeed at ml',       'wind_ml',   'm s-2',         'wind_speed') #(longname,shortname,units,stdname)
            self.names['winddir_ml'] = ('Wind direction at ml',  'winddir_ml','degrees north', 'wind_from_direction') #(longname,shortname,units,stdname)
        else:
            self.names['Vamp10m'] = ('Windspeed at 10m',      'Vamp10m',   'm s-2',         'wind_speed') #(longname,shortname,units,stdname)
            self.names['Vdir10m'] = ('Wind direction at 10m', 'Vdir10m','degrees north', 'wind_from_direction') #(longname,shortname,units,stdname)
            self.names['Vamp']    = ('Windspeed at ml',       'Vamp',      'm s-2',         'wind_speed') #(longname,shortname,units,stdname)
            self.names['Vdir']    = ('Wind direction at ml',  'Vdir',      'degrees north', 'wind_from_direction') #(longname,shortname,units,stdname)

        self.names['gust']       = ('Windgust at 10m',       'gust10',    'm s-2',         'wind_speed') #(longname,shortname,units,stdname)
        self.names['gustdir']    = ('Gust direction at 10m', 'gustdir10', 'degrees north', 'wind_from_direction') #(longname,shortname,units,stdname)


        for s in self.stationlist.keys():

            
            u10   = self.files[s].variables['u10m']
            v10   = self.files[s].variables['v10m']
            ugust = self.files[s].variables['u10mg']
            vgust = self.files[s].variables['v10mg']
            u_ml  = self.files[s].variables['u']
            v_ml  = self.files[s].variables['v']

            if 0: # naming
                self.create1D(self.files[s],'wind10',   towindspeed(u10[:],v10[:]))
                self.create1D(self.files[s],'winddir10',towinddir(u10[:],v10[:])  )
                self.create2D(self.files[s],'wind_ml',  towindspeed(u_ml[:],v_ml[:])  )
                self.create2D(self.files[s],'winddir_ml',towinddir(u_ml[:],v_ml[:])   )
            else:
                self.create1D(s,self.files[s],'Vamp10m', towindspeed(u10[0,:],v10[0,:]))
                self.create1D(s,self.files[s],'Vdir10m', towinddir(u10[0,:],v10[0,:])  )
                self.create2D(s,self.files[s],'Vamp',    towindspeed(u_ml[0,:],v_ml[0,:])  )
                self.create2D(s,self.files[s],'Vdir',    towinddir(u_ml[0,:],v_ml[0,:])   )

            self.create1D(s,self.files[s],'gust',     towindspeed(ugust[0,:],vgust[0,:]))
            self.create1D(s,self.files[s],'gustdir',  towinddir(  ugust[0,:],vgust[0,:]))


    def postproc_pz(self):

        # use a,b and surface pressure to calculate p at half levels, then full levels:
        a = self.ab[:self.nlevs+1 ]
        b = self.ab[ self.nlevs+1:]

        # give variables same dimension using numpy.tile:
        hlevels = self.nlevs + 1
        steps   = self.info['ltmax']+1

        self.names['p']     = ('Pressure at ml','p', 'hPa', 'air_pressure')
        self.names['height_f']     = ('Height at ml'  ,'height_f',   'm', 'height')

        # value taken from: McIlveen, Fundamentals of weather and Climate, p.74
        R_dry_air = 287.    # [J/(K.kg)]
        delta     = 1.e-10  # small number

        for s in self.stationlist.keys():

            #p_surf  = self.nc.variables['Psurf'][0,:]
            p_surf  = self.files[s].variables['Psurf'][0,:]

            ## calculate pressure at half, and full levels
            psml    = np.tile(p_surf[:],(hlevels,1)).T
            aml,bml = np.tile(a,(steps,1)),np.tile(b,(steps,1))
            p_half  = aml + psml * bml
            p_full  = (p_half[:,:-1] + p_half[:,1:])/2.

            phi_hl  = np.zeros((self.info['ltmax']+1,self.nlevs+1)) # allocate a new array
            ln_p_hl = np.log(p_half + delta) # calculate log of pressure

            # write as nc variable
            p    = self.create2D(s,self.files[s],'p',p_full)

            ## calculate thickness/height of levels: 
            t       = self.files[s].variables['t'][0,:]
            q       = self.files[s].variables['q'][0,:]
            tv      = t + (1 + 0.61 * q)   # virtual temperature
            try:
                oro     = self.files[s].variables['z_oro'][0,:]
            except:
                print 'No orography in netCDF, take station height'
                self.names['z_oro']     = ('Orography','z_oro', 'm','height')
                oro = np.array([stationlist[s]['z'] for step in range(steps)])
                self.create1D(s,self.files[s],'z_oro',oro)

        
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
            z_hl  = phi_hl/self.standard_g 
            z_mid = (z_hl[:,:-1] + z_hl[:,1:])/2
            self.create2D(s,self.files[s],'height_f', z_mid)


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

        def calc_td(t,q,p):
            tt   = t[:] - t_0
            tt[t == 0] = 1.e-10

            if 1: # test
                esat = 0.61078 * np.exp(17.2694 * (tt)/(tt + t_0 - 35.86)) # Bij Cisco 6.122 als factor?
                qsat = 0.622 * esat/((1.e-3)*p[:] - esat)
                rh   = q[:]/qsat 
                #td   = t_0 + (243.5*np.log(esat * rh/6.112))/(17.67- np.log(esat * rh/6.112))

                a,b,c = 6.112, 17.67, 243.5
                ps = a * np.exp( b * tt/ (c + tt))
                pa = rh * ps
                pa[pa<=0] = 1.e-10
                td = c * np.log(pa/a)/(b - np.log(pa/a))

            else: #except:
                print 'Rh, dew point calculation derailed'

            return rh,td
            
        self.names['td']     = ('Dew point temperature at ml','td', 'K', 'dew_point_temperature')
        self.names['rh']     = ('Relative humidity at ml',    'rh', '%', 'relative humidity')

        for s in self.stationlist.keys():

            t       = self.files[s].variables['t'][0,:]
            q       = self.files[s].variables['q'][0,:]
            p       = self.files[s].variables['p'][0,:] # created in postproc_pz

            rh,td   = calc_td(t,q,p)

            self.create2D(s,self.files[s],'td', td + t_0)
            self.create2D(s,self.files[s],'rh',rh)

    def postproc_percent(self):

        for s in self.stationlist.keys():
            for percvar in ['rh','rh2m','lcc','mcc','hcc','tcc','cc_a']:
                self.files[s].variables[percvar][:] = 100 * self.files[s].variables[percvar][:]

    def wrapup(self):
        
        for s in self.stationlist.keys():

            # test of dit uitmaakt:
            vstation = self.files[s].createVariable('stationid','u8',('station',)) 
            #vstation = self.files[s].createVariable('station','i4',('station',)) 
            vstation[0] = str(self.stationlist[s]['sid'])
            vstation.long_name = 'Station id'
            vstation.cf_role   = 'timeseries_id'


            self.files[s].close()



        print 'In nc_class: wrapped up'
        report_mem_usage() # report
