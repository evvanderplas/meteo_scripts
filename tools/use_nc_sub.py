#! /usr/bin/env python

import os,sys
import datetime as dt
import numpy as np
import subprocess
import cPickle as pickle

#import argparse, 
import pygrib
import netCDF4
import time
from netCDF4 import Dataset, date2num, num2date

from scipy.interpolate import griddata

from nc_profile_class import ncprofile, get_varname, stationlist

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

    from scipy.spatial import cKDTree
    
    mytree = cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)

    return indexes


def find_four_neighbours(latlons,(xlat,xlon),nearest=False,verb=True):

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
    if nearest: 
        if 0 not in (i0,j0): return (j0,i0)
        else: return None

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

def setup_nc(filename,startdate=None,station='Cabauw',stationid=6348,use_class = True):

    if use_class and startdate:
        netcdf = ncprofile(filename,startdate=startdate,station=station, stationid=stationid)

    else:
        ## set up a netCDF file
        netcdf = Dataset( filename, 'w', file_format='NETCDF4' ) 
        netcdf.set_fill_off()

        # create dims
        dtime = netcdf.createDimension('time', None)
        
        # create vars
        vtime = netcdf.createVariable('time','i4',('time',))
        #vtime.units = 'hours since 1950-01-01 00:00:00'
        vtime.units = 'seconds since 1970-01-01 00:00:00'

    return netcdf

def init_from_grib(ncObject,grb,use_class=True):

    nlevs = int(grb.numberOfVerticalCoordinateValues/2 -1)
    pv    = grb.pv
    lats,lons = grb.latlons()
    if use_class:
        ncObject.initpv(nlevs,pv,(lats,lons))
    else:

        dlevel = ncObject.createDimension('level', nlevs)
        hlevel = self.nc.createDimension('hlevel', self.nlevs+1)
        plevel = ncObject.createDimension('plevel', nlevs+1)

        a    = ncObject.createVariable('a','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
        a[:] = pv[:nlevs+1]

        b    = ncObject.createVariable('b','f4',('plevel',),zlib=True,fill_value=netCDF4.default_fillvals['f4'])   
        b[:] = pv[nlevs+1:]
    
    return ncObject


def init_index(grb,stationlist): #latlons):

    '''
    takes latlons from grib and the globally defined "stationlist" 
    to find indices of four surrounding points, to interpolate
    '''
    
    latlons = grb.latlons()

    # init indices, ONCE
    for s in stationlist.keys():
        print 'find index for ',stationlist[s]['name'] 
        
        xlat,xlon = stationlist[s]['latlon']
        stationlist[s]['ind'] = find_four_neighbours(latlons,(xlat,xlon),nearest=True,verb=True) #False)
        if stationlist[s]['ind'] == None: 
            del(stationlist[s])
            print stationlist 

    return latlons,stationlist # unnecessary? Does it update the global dict()?

def get_varname_local(parameter,level,leveltype,gribtab,model='harmonie'):

    if model == 'harmonie':
        field_codes = field_codes_harmonie


    if level_types.has_key(leveltype):
        ltype = level_types[leveltype]
    else: # default
        ltype = '105'   

    # if modellevel: search in index for level 001:
    if ltype == '109': 
        ilevel = 1
    # if pressure level: search for level 100 (?)
    elif ltype == '105' and int(level) in (100,250,300,500,800,850,900,925,950,1000):
        print 'pressure level? ',leveltype
        ilevel = 100
    else:
        ilevel = level

    
    # construct index:
    index= "{t:03}-{p:03}-{l:03}-{lt}".format(t=gribtab,p=parameter,l=ilevel,lt=ltype)
    if field_codes.has_key(index):
        (longname, shortname, units) = field_codes[ index ]
    else:
        # for now, just take vars that are described in table
        print 'Not in field_codes:',index
        return None,None,None

    if shortname is None or shortname.lower() == 'unknown':
        shortname = index
    elif ltype == '109': 
        shortname += '_ml'
                
    if longname is None:
        longname = grb.longName
    if units is None:
        units = '?'
    
    print 'From ',parameter,leveltype,level,' distilled ',longname,shortname,units
    return longname,shortname,units

def ncvar(ncObject,longname,shortname,units):

    #print ncObject.variables
    if shortname in ncObject.variables:
        vgrib = ncObject.variables[shortname]
        return vgrib

    elif shortname[-3:] == '_ml': # str(ltype) == '109':
        vgrib = ncObject.createVariable(shortname,'f4',('time','level'),
                                      zlib=True,
                                      fill_value=netCDF4.default_fillvals['f4'])

    #elif str(ltype) == '105' and grb.typeOfLevel == 'heightAboveGround': 
    elif 1: #grb.typeOfLevel in ('heightAboveGround','heightAboveSea'): 
        vgrib = ncObject.createVariable(shortname,'f4',('time',),
                                      zlib=True,
                                      fill_value=netCDF4.default_fillvals['f4']) 

    else:
        print 'Don t know yet:',longname,shortname
        
    vgrib.longname = longname
    vgrib.coordinates = "lat lon"
    #vgrib.location = stationlist[s]['name']
    vgrib.units    = units
    #vgrib.param    = par
    #vgrib.ltype    = ltype

    return vgrib
  
def add_var_to_nc(ncObject,grb,h,stationlist,latlons=None):

    '''
    take a netCDF handler, a grib message, leadtime h and a stationlist, (open grb only once)
    extract info and interpolated value, 
    add to netcdf file
    '''

    par        = grb.indicatorOfParameter
    leveltype  = grb.indicatorOfTypeOfLevel
    levelctype = grb.typeOfLevel
    level      = grb.level
    gribtab    = grb.table2Version

    longname,shortname,units = get_varname(par,level,levelctype,gribtab,model='harmonie')
    if longname == None: 
        print 'No value added to netCDF file ',par,levelctype,level 
        return ncObject

    vals    = grb.values
    try:
        latlons[0].shape == vals.shape
    except:
        latlons = grb.latlons()
           
    print type(stationlist), stationlist

    # actual interpolation
    for s in stationlist.keys():
        
        xlat,xlon = stationlist[s]['latlon']
        
        # do the interpolation
        #print 'interpolating for station ',s,stationlist[s]['name']
        iv = intValues(vals,latlons,stationlist[s]['ind'],(xlat,xlon),verb=False)
               
        # get or create/init the corresponding netCDF variable
        if type(ncObject == dict):
            vgrib = ncvar(ncObject[s],longname,shortname,units)
        else:
            print 'Waring: convention is 1 location per netcdf file (Feb 2014)'
            vgrib = ncvar(ncObject,longname,shortname,units)

        if shortname[-3:] == '_ml': #try: # if 2D field
            print 'writing ml: ',s,shortname,h,level-1
            vgrib[h,level-1] = iv
        elif shortname[-3:] == '_pl': 
            print 'writing pl: ',s,shortname,h,level-1
            vgrib[h,level-1] = iv
        else: #except:
            try: # if 1D field
                vgrib[h] = iv
            except: 
                print "Error with variable: ", longname #break
        if 0:#finally:
            #print 'done'
            pass

    # return the object for the next grib message
    return ncObject
    




stationlist_local = {
    1:{'name':'Cabauw',
       'latlon':(51.971, 4.927),
       },
    2:{'name':'Chilbolton',
       'latlon':(51.1445, -1.437),
       },
    3:{'name':'Palaiseau',
       'latlon':(48.713, 2.204),
       },
    4:{'name':'Lindenberg',
       'latlon':(52.1, 14.3),
       },
    5:{'name':'Juelich',
       'latlon':(50.909237, 6.413622),
       },
    6:{'name':'Schiphol',
       'latlon':(52.3156,4.7903), 
       },
    7:{'name':'De_Bilt',
       'latlon':(52.0989,5.1797), 
       },
    }

mlpars = np.array([   11, 33, 34, 40, 51, 58, 71, 76, 181, 184, 200, 201, 212, 213],dtype='int')
plpars = np.array([6, 11, 33, 34, 52],dtype='int') # 51 sometimes?
pressurelevels = [1000., 950., 925., 900., 850., 800., 700., 600., 500., 400., 300., 200.]


if __name__ == '__main__':

    print 'converting grib to netCDF4 profile using Python'

    gdir = '/nobackup_1/users/plas/temp'
    gdir = '/home/plas/data/vetrefo'
    fformat = 'HARM_N25_{sdtg}00_{lt}00_GB'
    field_codes = field_codes_harmonie


    sdate = dt.datetime(2012,10,25,0)
    sdate = dt.datetime(2010,7,14,6)
    #for h in hour_range(sdate,leadtime=48):

    subcmd  = 'nice -5 python subnc_slave.py {f} {sdtg} {lt}'
    pklfile = 'profarr_{sdtg}_{lt}.pkl'
    last    = 30
    filldict = {
        'model':'HARM_REFO',
        'sdtg':sdate.strftime('%Y%m%d%H'),
        'lt':str(last).zfill(3),
        }

    # get one gribmessage to initialise:
    gfile  = os.path.join(gdir,fformat.format(**filldict))
    grbs   = pygrib.open(gfile)
    grbref = grbs[1]
    grbs.close()

    # interpolate once (!)
    latlons,stationlist = init_index(grbref,stationlist)

    # dump info to communicate with threads:
    sf = open('statlist.pkl','wb')
    pickle.dump(stationlist,sf)
    sf.close()

    # setup a netCDF file for each station:
    for s in stationlist:
        fname  = 'profile_{model}_{loc}_{sdtg}.nc'.format(loc=stationlist[s]['name'],**filldict)
        stationlist[s]['nc'] = setup_nc(fname,startdate=sdate,
                                        station=stationlist[s]['name'],
                                        stationid=stationlist[s]['sid'])
        stationlist[s]['nc'] = init_from_grib(stationlist[s]['nc'],grbref)

    # if last pickle file does not (yet) exist: parse gribfiles:
    lastpkl = pklfile.format(**filldict)
    print 'Check: need to run?',lastpkl,os.path.exists(lastpkl)
    if not os.path.exists(pklfile.format(**filldict)):
        newrun=True
        pid = {}
        for h in range(0,last+1,1):
            filldict['lt'] = str(h).zfill(3)
            gfile = os.path.join(gdir,fformat.format(**filldict))
            print gfile,  os.path.exists(gfile)
            pid[h] = subprocess.Popen(subcmd.format(f=gfile,**filldict),shell=True)

            # test
            pid[h].communicate()
    else:
        newrun=False


    nroflevels   = 60
    nrofstations = max(stationlist.keys())

    mlprof = {}
    for p in mlpars:
        mlprof[p] = np.zeros((nroflevels,5+nrofstations,last),dtype='float')

    plprof = {}
    for p in plpars:
        plprof[p] = np.zeros((len(pressurelevels),5+nrofstations,last),dtype='float')

    surfprof = {}

    prof = {}
    for h in range(0,last,1):
        filldict['lt'] = str(h).zfill(3)

        if newrun: pid[h].communicate()
        
        # picklefile contans a numpy array 
        # (param,leveltype,level,timerange,value)
        prof[h] = pickle.load(open(pklfile.format(**filldict),"rb"))

        for p in mlpars:
            u = prof[h][(prof[h][:,0]==p) & (prof[h][:,1] == 109),:]
            #u.view('i4,i4,i4,i4,f4').sort(order=['f2'],axis=0)
            mlprof[p][:,:,h] = u 

        for p in plpars:
            u = prof[h][(prof[h][:,0]==p) & (prof[h][:,1] == 100),:]
            #u.view('i4,i4,i4,i4,f4').sort(order=['f2'],axis=0)
            plprof[p][:,:,h] = u 


        for item in prof[h][(prof[h][:,1] != 100) & (prof[h][:,1] != 109),:]:
            p,lt,l,tr,tab = item[0:5]
            index = "{t:03}-{p:03}-{l:03}-{lt}".format(t=int(tab),p=int(p),l=int(l),lt=int(lt))
            try:
                surfprof[index][:,h] = item
            except:
                print 'Make array first'
                surfprof[index] = np.zeros((5+nrofstations,last),dtype='float')
                surfprof[index][:,h] = item

        #u.view('i4,i4,i4,i4,f4').sort(order=['f2'],axis=0)
    print mlprof[11][:,7,:]


    for p in mlpars:
        gribtab = mlprof[p][0,4,0]
        print 'getting variable name',int(p),1,109,int(gribtab)
        longname,shortname,units,standardname = get_varname(int(p),1,'hybrid',int(gribtab),model='harmonie',version=3712)
        if longname == None: 
            print 'Handle exception'
            continue #sys.exit(0)
        for s in stationlist.keys():
            stationlist[s]['nc'].create2D(shortname,longname,units,standardname)
            stationlist[s]['nc'].nc.variables[shortname][:,:] = mlprof[p][:,5+s-1,:].T

    for p in plpars:
        gribtab = plprof[p][0,4,0]
        print 'getting variable name',int(p),1,100,int(gribtab)
        longname,shortname,units,standardname = get_varname(int(p),1,'isobaricInhPa',int(gribtab),model='harmonie',version=3712)
        if longname == None: 
            print 'Handle exception'
            sys.exit(0)
        for s in stationlist.keys():
            stationlist[s]['nc'].create_pl(shortname,longname,units,standardname)
            stationlist[s]['nc'].nc.variables[shortname][:,:] = plprof[p][:,5+s-1,:].T

        #for h in range(0,last,1):
    for p,leveltype,level,tr,tab in prof[h][(prof[h][:,1] != 100) & (prof[h][:,1] !=109),:][:,:5]:
        #print p,leveltype,level,tr,gribtab
        index = "{t:03}-{p:03}-{l:03}-{lt}".format(t=int(tab),p=int(p),l=int(level),lt=int(leveltype))
        longname,shortname,units,standardname = get_varname(int(p),int(level),int(leveltype),int(gribtab),model='harmonie',version=3712)
        if longname == None: 
            print 'Handle exception'
            continue #sys.exit(0)
        for s in stationlist.keys():
            if shortname not in stationlist[s]['nc'].nc.variables.keys():
                stationlist[s]['nc'].create1D(shortname,longname,units,standardname)
            stationlist[s]['nc'].nc.variables[shortname][:] = surfprof[index][5+s-1,:]
            if index == '253-011-002-105': print surfprof[index][5+s-1,:]
           
    sys.exit(0)

    indices = {}
    for h in range(0,4,1):
    #for h in range(0,49,1):
        gfile = os.path.join(gdir,fformat.format(dtgh = sdate.strftime('%Y%m%d%H'),lt = str(h).zfill(3)))
        print h,gfile, os.path.exists(gfile)

        grbs = pygrib.open(gfile)
        try:
            if init: pass
        except NameError:
            ncprof              = init_from_grib(ncprof,grbs[1])

            print type(stationlist), stationlist
            latlons,stationlist = init_index(grbs[1],stationlist)
            print type(stationlist), stationlist

            init = True

        # insert a time
        vtime = ncprof.variables['time']
        writetime = int( date2num(sdate + dt.timedelta(hours = h,seconds=30), 
                                  vtime.units, 
                                  calendar='standard'))
        vtime[h] = writetime

        # loop over fields
        for grb in grbs:

            add_var_to_nc(ncprof,grb,h,stationlist,latlons)

        # close the grib file
        grbs.close()

    # close the netCDF file, exit
    ncprof.close()
    sys.exit(0)
