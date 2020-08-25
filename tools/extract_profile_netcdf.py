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

from harmoniefields import level_types,field_codes_harmonie

def hour_range(begindate,leadtime = 48,dl=1,enddate = None):

    if type(begindate) != dt.datetime: 
        print 'hour_range should be called with datetime object, got',begindate
        sys.exit(1)

    begindate = begindate.replace(minute=0,second=0,microsecond=0)

    for h in range(0,leadtime,dl):
        yield begindate + dt.timedelta(hours=h)

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

def setup_nc(filename):

    ## set up a netCDF file
    netcdf = Dataset( filename, 'w', file_format='NETCDF4' ) 
    netcdf.set_fill_off()

    # create dims
    dtime = netcdf.createDimension('time', None)
    
    # create vars
    vtime = netcdf.createVariable('time','i4',('time',))
    vtime.units = 'hours since 1950-01-01 00:00:00'

    return netcdf

def init_from_grib(ncObject,grb):

    nlevs = int(grb.numberOfVerticalCoordinateValues/2 -1)
    dlevel = ncObject.createDimension('level', nlevs)
    plevel = ncObject.createDimension('plevel', nlevs+1)

    pv     = grb.pv
    a    = ncObject.createVariable('a','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
    a[:] = pv[:nlevs+1]

    b    = ncObject.createVariable('b','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
    b[:] = pv[nlevs+1:]
    
    return ncObject


def init_index(latlons):

    '''
    takes latlons and the globally defined "stationlist" 
    to find indices of four surrounding points, to interpolate
    '''
    
    # init indices, ONCE
    for s in stationlist.keys():
        print 'find index for ',stationlist[s]['name'] 
        
        xlat,xlon = stationlist[s]['latlon']
        stationlist[s]['ind'] = find_four_neighbours(latlons,(xlat,xlon),verb=False)

    return stationlist # unnecessary? Does it update the global dict()?

def get_varname(parameter,level,leveltype,gribtab,model='harmonie'):

    # if modellevel: search in index for level 001:
    if ltype == '109': 
        ilevel = 1
    else:
        ilevel = level

    
    # construct index:
    index= "{t:03}-{p:03}-{l:03}-{lt}".format(t=gribtab,p=par,l=ilevel,lt=ltype)
    if field_codes.has_key(index):
        (longname, shortname, units) = field_codes[ index ]
    else:
        # for now, just take vars that are described in table
        print 'Not in field_codes:',index
        return None,None,None

    if shortname is None or shortname.lower() == 'unknown':
        shortname = index
                
    if longname is None:
        longname = grb.longName
    if units is None:
        units = '?'
    
    print 'From ',parameter,leveltype,level,' distilled ',longname,shortname,units
    return longname,shortname,units

def add_var_to_nc(ncObject,grb):

    '''
    take a grib message, extract info and interpolated value, 
    add to netcdf file
    '''

    return ncObject
    


stationlist = {
    1:{'name':'Cabauw',
       'latlon':(51.57889, 4.53872),
       }
    }



if __name__ == '__main__':

    print 'converting profile to netCDF4 using Python'

    gdir = '/nobackup_1/users/plas/temp'
    fformat = 'HARM_N25_{dtgh}00_{lt}00_GB'
    field_codes = field_codes_harmonie

    ncprof = setup_nc('new_ncfile.nc')

    ## set up a netCDF file
    netcdf = Dataset( 'profile.nc', 'w', file_format='NETCDF4' ) 
    netcdf.set_fill_off()

    # create dims
    dtime = netcdf.createDimension('time', None)
    #dlevel= netcdf.createDimension('level', 60) # later
    #dlat  = netcdf.createDimension('rlat', 1)
    #dlon  = netcdf.createDimension('rlon', 1)

    # create vars
    #vlats = netcdf.createVariable('lat','f4',('rlat','rlon'))
    #vlons = netcdf.createVariable('lon','f4',('rlat','rlon'))
    #vrlons = netcdf.createVariable('rlon','f4',('rlon'))
    #vrlats = netcdf.createVariable('rlat','f4',('rlat'))
    vtime = netcdf.createVariable('time','i4',('time',))
    vtime.units = 'hours since 1950-01-01 00:00:00'
    #vproj = netcdf.createVariable('projection','c', () )


    sdate = dt.datetime(2012,10,25,0)
    #for h in hour_range(sdate,leadtime=48):


    indices = {}
    for h in range(0,4,1):
        gfile = os.path.join(gdir,fformat.format(dtgh = sdate.strftime('%Y%m%d%H'),lt = str(h).zfill(3)))
        print h,gfile, os.path.exists(gfile)

        grbs = pygrib.open(gfile)
        try:
            print 'filled properties ',nlevs
        except NameError:
            ncprof = init_from_grib(ncprof,grbs[1])

            nlevs = int(grbs[1].numberOfVerticalCoordinateValues/2 -1)
            dlevel= netcdf.createDimension('level', nlevs)
            plevel= netcdf.createDimension('plevel', nlevs+1)

            print 'nr of levs ',nlevs
            pv    = grbs[1].pv
            print pv
            print 'nr of as,bs ',len(pv[:nlevs+1]),len(pv[nlevs+1:])

            a = netcdf.createVariable('a','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
            a[:] = pv[:nlevs+1]

            b = netcdf.createVariable('b','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
            b[:] = pv[nlevs+1:]


        # insert a time
        vtime = ncprof.variables['time']
        writetime = int( date2num(sdate + dt.timedelta(hours = h,seconds=30), vtime.units, calendar='standard'))
        vtime[h] = writetime


        for grb in grbs:

            par = grb.indicatorOfParameter
            leveltype  = grb.indicatorOfTypeOfLevel
            levelctype = grb.typeOfLevel
            level = grb.level
            gtime = sdate + dt.timedelta(hours=h)

            if level_types.has_key(grb.typeOfLevel):
                ltype = level_types[grb.typeOfLevel]
            else: # default
                ltype = '105'
                
            #print grb.table2Version,par,level,ltype,type(ltype)

            # if modellevel: search in index for level 001:
            if ltype == '109': 
                ilevel = 1
            else:
                ilevel = level

            index= "{t:03}-{p:03}-{l:03}-{lt}".format(t=grb.table2Version,p= par, l=ilevel, lt=ltype)
            if field_codes.has_key(index):
                (longname, shortname, units) = field_codes[ index ]
            else:
                # for now, just take vars that are described in table
                print 'Not in field_codes:',index
                continue
            if shortname is None or shortname.lower() == 'unknown':
                shortname = index
                #if grb.shortName == 'unkown':
                #    shortname = index #str(par).zfill(3)
                #else:
                #    shortname = grb.shortName
                
            if longname is None:
                longname = grb.longName
            if units is None:
                units = '?'

            #print par,leveltype,levelctype,ltype,level,longname,shortname,gtime.strftime('%Y%m%d%H')
            #print 'units',units,type(units)

            vals = grb.values
            try:
                if vals.shape != latlons[0].shape:
                    print 'Wrong message: lats,lons do not coincide with values: ',vals.shape,latlons[0].shape
            except NameError:
                latlons = grb.latlons()

            # this
            stationlist = init_index(latlons)
            
            # for that
            # init indices, ONCE
            for s in stationlist.keys():

                try: 
                    #print 'Indices for ',stationlist[s]['name'],stationlist[s]['ind'] 
                    stationIndex = stationlist[s]['ind']
                except KeyError:
                    print 'find index for ',stationlist[s]['name'] 
                    
                    xlat,xlon = stationlist[s]['latlon']
                    stationlist[s]['ind'] = find_four_neighbours(latlons,(xlat,xlon),verb=False)


            # actual interpolation
            for s in stationlist.keys():
                
                xlat,xlon = stationlist[s]['latlon']

                #indices = find_four_neighbours(latlons,(xlat,xlon),verb=False)
                iv = intValues(vals,latlons,stationlist[s]['ind'],(xlat,xlon),verb=False)
                print grb.shortName, shortname, iv

                if ltype == '109':
                    #
                    # multilevel field on hybrid levels
                    #
                    #continue # EvdP focus on single
                    shortname = shortname+'_ml'
                    if shortname in netcdf.variables:
                        vgrib = netcdf.variables[shortname]
                    else:
                        vgrib = netcdf.createVariable(shortname,'f4',('time','level'),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
                        #vgrib = netcdf.createVariable(shortname,'f4',('level','time'),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
                        vgrib.longname = longname
                        vgrib.coordinates = "lat lon"
                        vgrib.location = stationlist[s]['name']
                        #vgrib.grid_mapping = "projection"
                        vgrib.units    = units
                        vgrib.param    = par
                        vgrib.ltype    = ltype

                    print 'writing ml: ',shortname,h,level-1
                    try:
                        vgrib[h,level-1] = iv
                    except ValueError:
                        print "Error with variable: ", longname #break

                elif str(ltype) == '105' and grb.typeOfLevel == 'heightAboveGround': 
                    #
                    # single level field
                    # 
                    #continue # EvdP focus on ml
                    print ltype, grb.typeOfLevel,grb

                    if shortname in netcdf.variables:
                        vgrib = netcdf.variables[shortname]
                    else:
                        vgrib = netcdf.createVariable(shortname,'f4',('time',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])
                        vgrib.longname = longname
                        vgrib.units    = units
                        vgrib.param    = par
                        vgrib.ltype    = ltype

                    if shortname == 'u': print 5*'=',h,'\n',grb


                    print 'writing single: ',shortname,h,writetime
                    try:
                        vgrib[h] = iv
                    except ValueError:
                        print "Error with variable: ", longname            #

        grbs.close()
    netcdf.close()
    sys.exit(0)

