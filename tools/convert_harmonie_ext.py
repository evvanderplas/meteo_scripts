#!/usr/bin/python

# Convert HARMONIE grib file to netCDF4 (compressed)
#
# Jisk Attema
# attema@knmi.nl
#
# last update: 10 June 2012 

import argparse, sys, os.path, re
import pygrib
import pyproj
import netCDF4
import time
from netCDF4 import Dataset #, date2num, num2date
import datetime as dt
#from datetime import datetime, timedelta
from numpy import *




# used later to check for missing values
maskarray = ma.zeros([1])


##################
#
# NB new harmonie codes!!
from harmonie_codes_kpt import *

def toposix(dtobj,starttime=dt.datetime(1970,1,1,0,0,0)):

    if type(dtobj) != dt.datetime:
        print 'something wrong with datetime:',dtobj
        sys.exit(1)

    timedl = dtobj - starttime

    return int(timedl.total_seconds())

def finddate(somefname):
    fc = re.compile('\S\D*(\d{4})(\d{2})(\d{2})(\d{2})\d*\D*(\d{3})\d*\D*$',re.VERBOSE)
    try:
        dtl = array(fc.search(somefname).groups(),dtype=int)
        #return dtl
        year,month,day,hour = dtl[0],dtl[1],dtl[2],dtl[3]
        leadtime = dtl[4]
        startdate = dt.datetime(year,month,day,hour,0,0) 
        validdate = dt.datetime(year,month,day,hour,0,0) + dt.timedelta(hours=int(leadtime))
        print 'found ',startdate,validdate
        #writetime = int( date2num( validdate, vtime.units, calendar='standard') )
        return startdate,validdate

    except:
        print 'Filename does not match yyyymmddhh bla lt bla',somefname
        return None,None


#
##################

##############################################################
#                Parse commandline
##############################################################

parser = argparse.ArgumentParser(description='Parse conversion arguments.')
parser.add_argument(      '--version',                              action='version', version='%(prog)s 1.0')


parser.add_argument('-O', '--overwrite',             default=False, action='store_true',
                    help="Overwrite existing file" )

parser.add_argument('-A', '--append',                default=False, action='store_true',
                    help="Append to existing file" )

parser.add_argument( '-n', '--nodouble',              default=False, action='store_true',
                    help="Do not overwrite fields even when grib param/level/type/tri are identical (conflicts with time series!)" )

parser.add_argument('-t', '--table',                 default='harmonie', action='store',
                    help='GRIB table (harmonie36 / harmonie37 / racmo)')

parser.add_argument('-s', '--nomult',                default=False, action='store_true',
                    help="Skip multilevel fields (surface only)" )

parser.add_argument('-a', '--nosurf',                default=False, action='store_true',
                    help="Skip single level fields (atmosphere only)" )

parser.add_argument('-o', '--output',                               action='store',      dest='outfile', required=True,
                    help="The netCDF4 ouput file" )

parser.add_argument('-l', '--levels',    type=int,   default=60,    action='store',      dest='nlevs',
                    help="The number of vertical levels" )

parser.add_argument('-i', '--index',                 default=False, action='store_true', dest='always_index',
                    help="Dont translate table/param/level/ltype." );

#parser.add_argument('-r', '--rtime',     type=str,   default='hours since 1950-01-01 00:00:00', action='store',
#                    help="NetCDF reference time." );
parser.add_argument('-r', '--rtime',     type=str,   default='seconds since 1970-01-01 00:00:00', action='store',
                    help="NetCDF reference time." );

parser.add_argument('-d', '--deltatime', type=int,   
                    #default=1,
                    default=3600,
                    help="Input file frequency in units of rtime." );

parser.add_argument('-v', '--verbose',               default=False, action='store_true',
                    help="Give more verbose output." );

parser.add_argument('-prt', '--date_from_filename',  default=False, action='store_true',
                    help="Infer datetime from filename. (try)" );

parser.add_argument('-m', '--metadata',  metavar='key=value',       action='append',     nargs='*',
                    help="key=value will be added to netCDF4 metadata" );

parser.add_argument('-k', '--keys',      metavar='TBL-PRM-LVL-LTP|shortname', action='append',     nargs='*',
                    help="Only process GRIB messages with given key. (for multi-level fields set LVL to 001)" );


parser.add_argument('infiles', metavar='inputfiles', type=str,      action='append',     nargs='*',
                    help='GRIB inputs file to convert')

args = parser.parse_args()

if args.verbose:
    print "Input:              ", args.infiles
    print "Table:              ", args.table
    print "Output:             ", args.outfile
    print "Overwrite:          ", args.overwrite
    print "Skip multilevel     ", args.nomult
    print "Skip single level   ", args.nosurf
    print "Levels:             ", args.nlevs
    print "Reference time:     ", args.rtime
    print "Name keys by index: ", args.always_index
    if args.keys:
        print "Keys:               ", args.keys
    else:
        print "Keys:                All"
    print "Metadata:           ",  args.metadata
    print "From filename:      ",  args.date_from_filename

#sys.exit(0)

# select table
write_fanames = False
if args.table == 'harmonie' or args.table == 'harmonie36':
    field_codes = field_codes_harmonie36
elif args.table == 'harmonie37':
    field_codes = field_codes_harmonie37
    write_fanames = True
elif args.table == 'racmo':
    field_codes = field_codes_racmo
else:
    assert False, "Unknown GRIB-table %s" % args.table

# check outputfile
if not args.outfile:
    assert False, "no output given"

if os.path.isfile( args.outfile ) and not args.overwrite and not args.append:
    assert False, "file exists, use -O to overwrite, or -A to append."

# check input file
if not args.infiles:
    assert False, "no input files"

flattend = []
if args.infiles:
    for sublist in args.infiles:
        flattend += sublist 
    args.infiles = flattend

# Check if we want a subset of keys
flattend = []
if args.keys:
    for extra in args.keys:
        flattend += extra

message_subset = {}
if flattend:
    for k in flattend:
        key_added=False
        if field_codes.has_key(k):
            message_subset[k] = True
            print "Adding key ", k
            continue
        else:
            # find the index for this shortname
            for sb in field_codes.keys():
                if field_codes[sb][1] == k:
                    message_subset[sb] = True
                    print "Adding key ", sb
                    key_added=True
        if not key_added:
            print "Could not find key ", k

##############################################################
#                Setup netcdf file
##############################################################


if args.append:
    # open a NetCDF file for output, dont prefill with missing values
    # missing values are handled later
    netcdf = Dataset( args.outfile, 'r+', file_format='NETCDF4' ) 
    netcdf.set_fill_off()
    dtime = netcdf.dimensions['time']
    vtime = netcdf.variables['time']
else:
    # open a NetCDF file for output, dont prefill with missing values
    # missing values are handled later
    netcdf = Dataset( args.outfile, 'w', file_format='NETCDF4' ) 
    netcdf.set_fill_off()

    # Load grid from second message form first file (RACMO-Asimov has a hidden first message, ignore it)
    grbs = pygrib.open(args.infiles[0])
    try:
        grb = grbs[2]
    except IOError:
        grb = grbs[1]

    griblats, griblons = grb.latlons()
    grbs.close()
    nlats, nlons = shape(griblats)

    ##
    # time (record) dimension
    #

    dtime = netcdf.createDimension('time', None)
    vtime = netcdf.createVariable('time','i4',('time',))

    vtime.units = args.rtime
    vtime.calendar = 'standard'
    vtime.long_name = 'time'
    vtime.standard_name = 'time'

    rtime = netcdf.createVariable('forecast_reference_time','i4',())
    rtime.units = args.rtime
    rtime.calendar = 'standard'
    rtime.long_name     = 'forecast_reference_time'
    rtime.standard_name = 'forecast_reference_time'

    startdate,validdate = finddate(args.infiles[0])
    #rtime[0] = int( date2num( startdate, rtime.units, calendar=rtime.calendar) )
    rtime[0] = toposix(startdate)
    

    ##
    # vertical dimension
    #
    if not args.nomult:
        dmlev = netcdf.createDimension('mlev', args.nlevs)
        dhyam = netcdf.createDimension('nhym', args.nlevs)
        dhyai = netcdf.createDimension('nhyi', args.nlevs + 1)

        vmlev = netcdf.createVariable('mlev','i4',('mlev'))
        vhyai = netcdf.createVariable('hyai','f4',('nhyi'))
        vhybi = netcdf.createVariable('hybi','f4',('nhyi'))
        vhyam = netcdf.createVariable('hyam','f4',('nhym'))
        vhybm = netcdf.createVariable('hybm','f4',('nhym'))

        # for harmonie, set ahalf and bhalf pressure level coordinates
        if re.match('^harmonie', args.table):
            if args.nlevs == 60:
                vhyai[:] = vertical_levels_harmonie_mf60[ 'AHALF' ]
                vhybi[:] = vertical_levels_harmonie_mf60[ 'BHALF' ]
                vhyam[:] = 0.5 * (vhyai[0:60] + vhyai[1:61] )
                vhybm[:] = 0.5 * (vhybi[0:60] + vhybi[1:61] )

        vmlev[:]  = range(args.nlevs)

        vmlev.standard_name = "hybrid_sigma_pressure" ;
        vmlev.long_name = "hybrid level at layer midpoints" ;
        vmlev.units = "level" ;
        vmlev.positive = "down" ;
        vmlev.formula = "hyam hybm (mlev=hyam+hybm*aps)" ;
        vmlev.formula_terms = "ap: hyam b: hybm ps: aps" ;

        vhyai.long_name = "hybrid A coefficient at layer interfaces" ;
        vhyai.units = "Pa" ;
        vhybi.long_name = "hybrid B coefficient at layer interfaces" ;
        vhybi.units = "1" ;
        vhyam.long_name = "hybrid A coefficient at layer midpoints" ;
        vhyam.units = "Pa" ;
        vhybm.long_name = "hybrid B coefficient at layer midpoints" ;
        vhybm.units = "1" ;

    ##
    # Tile dimension (RACMO)
    #

    if args.table == 'racmo':
        dtile = netcdf.createDimension('tile', 8 )

    ##
    # Horizontal grid
    # 

    print('Grid type: {}'.format(grb.gridType))

    # simple 2d regular lat lon grid
    #
    if grb.gridType in ['regular_ll','regular_gg']:

        model_grid_is_simple = True

        # Model grid, code taken from the pygrib sourceocde
        nx = grb.Ni
        ny = grb.Nj

        lon1 = grb.longitudeOfFirstGridPointInDegrees
        lon2 = grb.longitudeOfLastGridPointInDegrees

        if lon1 >= 0 and lon2 < 0 and grb.iDirectionIncrement > 0:
            lon2 = 360+lon2
        if lon1 >= 0 and lon2 < lon1 and grb.iDirectionIncrement > 0:
            lon1 = lon1-360

        lat1 = grb.latitudeOfFirstGridPointInDegrees
        lat2 = grb.latitudeOfLastGridPointInDegrees

        # workaround for grib_api bug with complex packing.
        # (distinctLongitudes, distinctLatitudes throws error,
        # so use np.linspace to define values)
        if grb.packingType.startswith('grid_complex'):
            # this is not strictly correct for gaussian grids,
            # but the error is very small.
            if lat1 < lat2:
                lats = np.linspace(lat1,lat2,ny)
            else:
                lats = np.linspace(lat2,lat1,ny)
                lons = np.linspace(lon1,lon2,nx)
        else:
            lats = grb.distinctLatitudes
            lons = grb.distinctLongitudes

        # don't trust distinctLongitudes
        # when longitudeOfLastGridPointInDegrees < 0
        # (bug in grib_api 1.9.16)
        if lon2 < 0:
            lons = np.linspace(lon1,lon2,nx)

        # Model grid is actual grid
        dlat  = netcdf.createDimension('lat', nlats)
        dlon  = netcdf.createDimension('lon', nlons)
        modelgrid = ('lat','lon', )

        vlons = netcdf.createVariable('lon','f4',('lon'))
        vlats = netcdf.createVariable('lat','f4',('lat'))

        vlats.units = 'degrees_north'
        vlats.long_name = 'latitude'
        vlats.standard_name = 'latitude'

        vlons.units = 'degrees_east'
        vlons.long_name = 'longitude'
        vlons.standard_name = 'longitude'

        vlats[:] = lats
        vlons[:] = lons

    else:

    # Complex grid, find model grid, projection, and actual latlons
    #
        model_grid_is_simple = False

        if grb.gridType in ('lambert',): #'rotated_ll'):

            # Model grid, code taken from the pygrib sourceocde
            lat1 = grb.latitudeOfFirstGridPointInDegrees
            lon1 = grb.longitudeOfFirstGridPointInDegrees
            try:
                nx = grb.Nx
                ny = grb.Ny
            except:
                nx = grb.Ni
                ny = grb.Nj
            dx = grb.DxInMetres
            dy = grb.DyInMetres
            pj = pyproj.Proj(grb.projparams)
            llcrnrx, llcrnry = pj(lon1,lat1)

            dlat  = netcdf.createDimension('y', nlats)
            dlon  = netcdf.createDimension('x', nlons)
            vrlats = netcdf.createVariable('y','f4',('y'))
            vrlons = netcdf.createVariable('x','f4',('x'))
            modelgrid = ('y','x',)

            vrlats.units = 'm'
            vrlats.long_name = 'y coordinate of projection'
            vrlats.standard_name = 'projection_y_coordinate'

            vrlons.units = 'm'
            vrlons.long_name = 'x coordinate of projection'
            vrlons.standard_name = 'projection_x_coordinate'

            vrlons[:] = llcrnrx+dx*arange(nx)
            vrlats[:] = llcrnry+dy*arange(ny)

            # projection description
            vprojname = 'Lambert_Conformal'
            vproj = netcdf.createVariable(vprojname,'i4', () )
            vproj.grid_mapping_name = "lambert_conformal_conic"
            vproj.standard_parallel = repr(grb.projparams["lat_1"]) + " " + repr(grb.projparams["lat_2"])
            vproj.longitude_of_central_meridian = repr(grb.projparams["lon_0"])
            vproj.latitude_of_projection_origin = repr(grb.projparams["lat_0"])


        elif grb.gridType in ['rotated_ll','rotated_gg']:

            # Model grid, code taken from the pygrib sourceocde
            # Model grid is actual grid
            dlat  = netcdf.createDimension('rlat', nlats)
            dlon  = netcdf.createDimension('rlon', nlons)

            vrlons = netcdf.createVariable('rlon','f4',('rlon'))
            vrlats = netcdf.createVariable('rlat','f4',('rlat'))
            modelgrid = ('rlat','rlon',)

            vrlats.units = 'degrees'
            vrlats.long_name = 'latitude in rotated pole grid'
            vrlats.standard_name = 'grid_latitude'

            vrlons.units = 'degrees'
            vrlons.long_name = 'longitude in rotated pole grid'
            vrlons.standard_name = 'grid_longitude'

            vrlons[:] = grb.distinctLongitudes
            vrlats[:] = grb.distinctLatitudes

            # projection description
            vprojname = 'rotated_pole'
            vproj = netcdf.createVariable(vprojname,'c', () )
            vproj.grid_mapping_name = "rotated_latitude_longitude"
            vproj.grid_north_pole_latitude = repr(grb.projparams["o_lat_p"])
            vproj.grid_north_pole_longitude = repr(grb.projparams["lon_0"] + 180.0)


        # Actual grid
        vlats  = netcdf.createVariable('lat','f4',modelgrid)
        vlons  = netcdf.createVariable('lon','f4',modelgrid)

        vlats.units = 'degrees_north'
        vlats.long_name = 'latitude'
        vlats.standard_name = 'latitude'
        vlats.grid_mapping = 'projection'

        vlons.units = 'degrees_east'
        vlons.long_name = 'longitude'
        vlons.standard_name = 'longitude'
        vlons.grid_mapping = 'projection'

        vlats[:] = griblats
        vlons[:] = griblons

        # Projection string
        projstring = "+proj=" + grb.projparams['proj']
        for key, value in grb.projparams.iteritems():
            if key != 'proj':
                if type(value) == type("string"):
                    projstring += " +" + key + "=" + value
                else:
                    projstring += " +" + key + "=" + repr(value)

        vproj.proj = projstring.strip()


# Bug in netcdf4: extending history doesnt work because we're stuck in datamode..
# update history
if 'history' in netcdf.ncattrs():
    history = netcdf.history
else:
    history = ""
netcdf.history = time.ctime(time.time()) + ": "+ ' '.join(sys.argv) + history 

# set all command-line provided metadata
if args.metadata:
    flattend = []
    for extra in args.metadata:
        flattend += extra

    for extra in flattend:
        try:
            key, value = extra.split('=')
        except ValueError: 
            key = extra
            value = ''
        exec "netcdf.%s = value" % key 


##############################################################
#                Convert messages
##############################################################

for infile in args.infiles:

    # set a time in netCDF file from the second grb message in the file (to skip racmo/asimov header)
    try:
        grbs = pygrib.open(infile)
    except IOError:
        print 'File already gone?!?&^*%'
        continue

    try:
        grb = grbs[2]
    except IOError:
        grb = grbs[1]

    if args.date_from_filename:
        fpath,fname = os.path.split(infile)
        startdate,validdate = finddate(fname)
        print fname,' so ',startdate,validdate
    else:
        year  = grb.validityDate / (100 * 100)
        month = (grb.validityDate - year * 100 * 100) / 100
        day   = grb.validityDate - year * 100 * 100 - month * 100
        hour  = grb.validityTime  / 100
        minute= grb.validityTime - hour * 100
        seconds= 30

        sdtg  = dt.datetime.strptime(str(grb.dataDate),'%Y%m%d') # yyyymmdd
        shour = grb.dataTime/100
        startdate = sdtg + dt.timedelta(hours=shour)

        print 20*'='
        print infile
        print 'date',grb.validityDate
        print 'time',grb.validityTime
        print list(vtime[...])
        print 'Check dtime, time',len(dtime),dt.datetime( year, month, day, hour, minute, seconds ); #continue
        validdate =  dt.datetime( year, month, day, hour, minute, seconds )
        print 'Check resulting datetime',

    #reftime   = int( date2num(startdate, vtime.units, calendar='standard') )
    #writetime = int( date2num(validdate, vtime.units, calendar='standard') )
    reftime   = toposix(startdate)
    writetime = toposix(validdate)

    #print reftime,writetime; sys.exit(0)

    # find the right index; the try/except/if/else are needed because netCDF/python throw a lot of exceptions and errors
    try:
        if len(dtime) > 0:
            if writetime is not None:
                windex = list(vtime[...]).index( writetime )
            else:
                windex = windex + 1
            #print 'Index for time if no ValueError:',windex
        else:
            windex = 0
        vtime[windex] = writetime
    except ValueError:
        
        if len(dtime) > 0:
            if writetime is not None:
                windex = (writetime - vtime[0]) / args.deltatime
            else:
                windex = windex + 1
            #print 'Index for time if ValueError:',windex
        else:
            windex = 0
        vtime[windex] = writetime
    print 'Index for time at last:',windex,vtime[:]; #grbs.close(); continue


    if args.verbose:
        #print "Parsing file '%s' for date %s, index %4i, value %i" % (infile, dt.datetime( year, month, day, hour, minute ), windex, writetime )
        print "Parsing file '%s' for index %4i, value %i" % (infile, windex, writetime )

    for grb in grbs:

        ssuffix = None
        sprefix = None

        lsuffix = None
        lprefix = None

        ##
        #
        # Format codes as ???
        #
        table = "%03i" % grb.table2Version
        param = "%03i" % grb.indicatorOfParameter
        tri  = grb.timeRangeIndicator
        if level_types.has_key(grb.typeOfLevel):
            ltype = level_types[grb.typeOfLevel]
        else:
            print grb.typeOfLevel
            #ltype = "%03i" % grb.typeOfLevel
            next

        if ltype == '109':

            # hybrid levels are listed under level 001 
            level = "001"


        #######################
        # manual bugfixes
        #######################

        elif args.table == 'racmo' and ltype == '001':

            # for racmo, tiled fields are listed under level 001
            level = "001"

        elif args.table == 'racmo' and ltype == '100':

            # fields at given hPa level
            level = "000"
            suffix = "%03i" % grb.level

        # for some leveltypes, GRIB allocates a different number of octets, leading to byteswaps or truncations

        elif (args.table == 'racmo'    and ltype == '102') :

            level = grb.level
            level = ((level & 0xFF) << 8) + ((level & 0xFF00) >> 8)
            level = "%03i" % level

        elif ( re.match( '^harmonie', args.table)  and ltype == '008') :

            if ( grb.level == 39424 ):
                level = 666
            elif ( grb.level == 39680 ):
                level = 667
            else:
                level = grb.level
            level = "%03i" % level

        else:
            level = "%03i" % grb.level

        ##############################################
        # Parse time length indicator for Harmonie37
        ##############################################
    
        if args.table == 'harmonie37':
            tri = grb.timeRangeIndicator
            if tri == 0:
                lprefix = 'Instantaneous '
            elif tri == 2:
                lsuffix = ' over interval'
            elif tri == 4:
                sprefix = 'a'
                lprefix = 'Accumulated '
            else:
                print "Cannot parse tri ", tri

        ##
        #
        # look up field-codes
        #
        index= "%s-%s-%s-%s" % (table, param, level, ltype)


        ##
        #
        # Verbose output
        #
        if args.verbose:
            if field_codes.has_key(index):
                print index, "(",grb.level,")", field_codes[ index ]
            else:
                print index

        ##
        #
        # skip multilevel fields?
        # EvdP include lowest modellevel, exclude pressure levels
        #

        #if args.nomult and ltype == '109' and level != args.nlevs:
        if args.nomult and ltype in ('100','109') and grb.level != args.nlevs:
            # ... yes.
            #print 'debug ',param,ltype,level,grb.level,args.nlevs
            continue

        ##
        #
        # skip single level fields?
        #
        if args.nosurf and ltype != '109':
            # ... yes.
            continue

        ##
        #
        # skip this message?
        #
        if message_subset and not message_subset.has_key(index):
            # ... yes.
            continue

        if not args.always_index and field_codes.has_key(index):

            longname  = field_codes[ index ][0]
            shortname = field_codes[ index ][1]
            units     = field_codes[ index ][2]

            if write_fanames:
                faname = field_codes[ index ][3]

            standard_name  = field_codes[ index ][4]


        else:
            longname = index
            shortname = index
            faname = "?"
            units = "?"
            standard_name = index

        #
        if not longname:
            longname = index
        if not shortname:
            shortname = longname
        if not shortname:
            shortname = index
        if not units:
            units = "?"
        if not standard_name:
            standard_name = index


        ##
        #
        # apply suffix / prefix
        #
        if sprefix:
            shortname = "%s%s" % (sprefix, shortname)
        if lprefix:
            longname = "%s%s"  % (lprefix, longname)

        if lsuffix:
            longname = "%s%s"  % (longname, lsuffix)
        if ssuffix:
            shortname = "%s%s" % (shortname, ssuffix)


        ##
        #
        # read values
        #
        values = grb.values
        if type(values) == type(maskarray):
            values = values.filled( netCDF4.default_fillvals['f4'] )

        ##
        #
        # Create variable and write
        # EvdP added possibility to add lowest model levels even if nomult = True
        # 

        #if ltype == '109':
        if ltype == '109' and not args.nomult:
            #
            # mulitlevel field on hybrid levels
            #
            if shortname in netcdf.variables:
                vgrib = netcdf.variables[shortname]
            else:
               vgrib = netcdf.createVariable(shortname,'f4',('time','mlev',) + modelgrid,zlib=True,fill_value=netCDF4.default_fillvals['f4'])
               vgrib.long_name = longname
               vgrib.standard_name = standard_name
               if not model_grid_is_simple:
                   vgrib.coordinates = "lat lon"
                   vgrib.grid_mapping = vprojname
               vgrib.units    = units
               vgrib.param    = param
               vgrib.ltype    = ltype
               vgrib.timeRange = grb.timeRangeIndicator
               if write_fanames:
                   vgrib.faname = faname 
            try:
                vgrib[windex,grb.level-1,...] = values
            except ValueError:
                print "Error with variable: ", longname, " ltype=",ltype

        elif args.table == 'racmo' and ltype == '001':
            #
            # ramco tiled fields
            #
            if shortname in netcdf.variables:
                vgrib = netcdf.variables[shortname]
            else:
               vgrib = netcdf.createVariable(shortname,'f4',('time','tile',) + modelgrid,zlib=True,fill_value=netCDF4.default_fillvals['f4'])
               vgrib.long_name = longname
               if not model_grid_is_simple:
                   vgrib.coordinates = "lat lon"
                   vgrib.grid_mapping = vprojname
               vgrib.units    = units
               vgrib.param    = param
               vgrib.ltype    = ltype
               vgrib.timeRange = grb.timeRangeIndicator
               if write_fanames:
                   vgrib.faname = faname 

            # netCDF starts at 0, tiles at 1

            if sys.byteorder == 'little':
                # tiles are stored big-endian but this cpu is little endian
                tile = grb.level
                tile = ((tile & 0xFF) << 8) + ((tile & 0xFF00) >> 8) - 1
            else:
                tile = grb.level - 1

            try:
                vgrib[windex,tile,...] = values
            except ValueError:
                print "Error with variable: ", longname, " ltype=",ltype
        else:
            #
            # single level
            #
            if shortname in netcdf.variables and args.nodouble:
               vgrib = netcdf.variables[shortname]
            else:
               if shortname in netcdf.variables:
                   print "Double message encountered\n"
                   while shortname in netcdf.variables:
                       shortname += "X"
                   print "Renaming to ", shortname
               vgrib = netcdf.createVariable(shortname,'f4',('time',) + modelgrid,zlib=True,fill_value=netCDF4.default_fillvals['f4'])
               vgrib.long_name = longname
               vgrib.standard_name = standard_name
               if not model_grid_is_simple:
                   vgrib.coordinates = "lat lon"
                   vgrib.grid_mapping = vprojname
               vgrib.units    = units
               vgrib.param    = param
               vgrib.ltype    = ltype
               vgrib.level    = level
               vgrib.timeRange = grb.timeRangeIndicator
               if write_fanames:
                   vgrib.faname = faname 

            try:
                vgrib[windex,...] = values
            except ValueError:
                print "Error with variable: ", longname, " ltype=", ltype


##############################################################
#                Clean up
##############################################################

netcdf.close()
grbs.close()
