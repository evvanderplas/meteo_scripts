#! /usr/bin/env python

import os,sys
import datetime as dt
import pygrib
from netCDF4 import Dataset, date2num, num2date

import use_ncprof as nc_extract
from use_ncprof import stationlist

#sys.path.append('/usr/people/plas/python/fabriek')

BULL_OPER_SA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULLSA',
    }

BULL_OPER_LA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL, large area)',
    'shortname': 'BULL_LA'
    }

sources = [BULL_OPER_SA]
#sources = [BULL_OPER_LA]

if __name__ == '__main__':

    ## create datetime of latest possible run:
    today = dt.datetime.today()
    today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0,microsecond=0)
    td = dt.datetime.strftime(today,'%Y%m%d%H')

    for src in sources:

        gdir   = src['dir']
        shname = src['shortname']
        
        for h in range(0,24,3):
            st = today - dt.timedelta(hours = h)
            
            print 5*'*','retry:',st

            ## set up a netCDF file for every location this run:
            ncprof,init = {},False
            for s in stationlist:
                location = stationlist[s]['name']
                outfile = 'profile_{m}_{ymdh}_{loc}.nc'.format(m=shname,ymdh = st.strftime('%Y%m%d%H'),loc = location)
                if os.path.exists(outfile): 
                    print 'huh?'; 
                    print 'delete! remake', outfile 
                    os.remove(outfile) # sys.exit(0)
                    ncprof[s] = nc_extract.setup_nc(outfile)

                else: #continue
                    ncprof[s] = nc_extract.setup_nc(outfile)
    
            #for lt in range(0,49,1):
            for lt in range(0,1,1):
                #sdate = st + dt.timedelta(hours=lt)
                gfile = os.path.join(gdir,src['fformat'].format(dtgh = st.strftime('%Y%m%d%H'),lt = str(lt).zfill(3)))
                print lt,gfile, os.path.exists(gfile)
                if not os.path.exists(gfile): break ## assume if one is missing, all are missing
                
                print 10*'=','Open ',gfile
                grbs = pygrib.open(gfile)
                if init: #try:
                    #if latlons and init: 
                    print 'Already got a ',ncprof,latlons[0].shape
                    pass
                else: #except NameError:
                    print 'Generating new netCDF dataset'
                    for s in stationlist: 
                        print stationlist[s]['name']
                        ncprof[s] = nc_extract.init_from_grib(ncprof[s],grbs[1])
                    latlons,stationlist     = nc_extract.init_index(grbs[1],stationlist)
                    init = True
                #finally:
                #    print 'initialisation failed,',sys.exit(0)
                    
                # insert a time
                for s in stationlist:
                    #print 'writing time in ',s,ncprof[s]
                    vtime = ncprof[s].variables['time']
                    writetime = int( date2num(st + dt.timedelta(hours = lt,seconds=30), 
                                              vtime.units, 
                                              calendar='standard'))
                    print 'writing time in ',s,writetime, ncprof[s],writetime                    
                    vtime[lt] = writetime

                for grb in grbs:
                    print grb
                    nc_extract.add_var_to_nc(ncprof,grb,lt,stationlist,latlons)
                    
                # close the grib file
                grbs.close()

            # close the netCDF file(s)
            for s in stationlist:
                print 'Closing',st,s,stationlist[s]['name']
                try:
                    ncprof[s].close()
                except KeyError:
                    print 'Not opened: not closed,',outfile

            # do a new init for the next iteration
            init = False
