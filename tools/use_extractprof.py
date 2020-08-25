#! /usr/bin/env python

import os,sys
import pygrib
import datetime as dt
from os.path import join as pjoin


import nc_profile_class as ncprof

BULL_OPER_SA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULLSA',
    'ltmax':48,
    }

BULL_OPER_LA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULLLA',
    'ltmax':48,
    }

VETREFO = {
    'dir' : '/nobackup_1/users/plas/temp/2010/07/14/06',
    #'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETREFO',
    'ltmax':48,
    }
VETSTD = {
    'dir' : '/nobackup_1/users/plas/temp/VETSTD',
    #'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETSTD',
    'ltmax':18,
    }
VETP2 = {
    'dir' : '/nobackup_1/users/plas/temp/VETP2',
    #'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETP2',
    'ltmax':18,
    }

stationlist = ncprof.stationlist
suites = (BULL_OPER_SA,)
suites = (VETREFO,)
suites = (VETSTD,VETP2)
suites = (VETP2,)
outdir = '/nobackup_1/users/plas/ncdat/prof'


if __name__ == '__main__':

    ## create datetime of latest possible run:
    today = dt.datetime.today()
    today = dt.datetime(2010,7,14,6)

    today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0,microsecond=0)
    td = dt.datetime.strftime(today,'%Y%m%d%H')

    for src in suites:

        gdir   = src['dir']
        shname = src['shortname']
        model,version = src['model'],src['version']
        outdir = pjoin(outdir,shname)
        try: 
            os.makedirs(outdir)
        except:
            print 'Output directory in place: ',outdir

        for h in range(0,24,3):
            st = today - dt.timedelta(hours = h)
            dtgh = st.strftime('%Y%m%d%H')
            ltmax = src['ltmax']
            fformat = src['fformat']
            lastfile = pjoin(gdir,fformat.format(dtgh = dtgh,lt = str(ltmax).zfill(3)))
            print lastfile, os.path.exists(lastfile)
            if not os.path.exists(lastfile): continue
            print 'so..'

            nlevs,ab,latlons = ncprof.get_levs_ab_latlon_from_grib(lastfile)
            stationlist      = ncprof.init_index(latlons,stationlist)

            ncperdate = {}
            for s in stationlist:
                location = stationlist[s]['name']
                outfile = 'profile_{m}_{dtgh}_{loc}.nc'.format(m=shname,dtgh = dtgh,loc = location)           
                #outfile = pjoin(outdir,shname,'{n}_')
                if os.path.exists(outfile): continue
                
                ncperdate[s] = ncprof.ncprofile(st,outfile,stationlist[s],stationlist[s]['sid'])
                ncperdate[s].initpv(nlevs,ab,latlons)
                
            for lt in range(0,ltmax+1,1):
            #for lt in range(0,2,1):

                # insert a time
                for  s in stationlist.keys():
                    ncperdate[s].writetime(st,lt)

                gfile = pjoin(gdir,fformat.format(dtgh = dtgh,lt = str(lt).zfill(3)))
                grbs = pygrib.open(gfile)
                for grb in grbs:
                    
                    par        = grb.indicatorOfParameter
                    leveltype  = grb.indicatorOfTypeOfLevel
                    levelctype = grb.typeOfLevel
                    level      = grb.level
                    gribtab    = grb.table2Version
                    
                    vals       = grb.values
  
                    print par,level,levelctype,gribtab,model,version

                    longname,shortname,units,stdname = ncprof.get_varname(par,
                                                                          level,
                                                                          levelctype,
                                                                          gribtab,
                                                                          model=model,
                                                                          version=version)

                    if shortname == None:
                        print 'Again: not in field codes'
                        continue

                    for s in stationlist.keys():
                        #xlat,xlon = ncperdate[s].stationlist[s]['latlon']
                        xlat,xlon = stationlist[s]['latlon']
                        iv = ncprof.intValues(vals,latlons,stationlist[s]['ind'],(xlat,xlon),verb=False)

                        ncperdate[s].add_value_to_var(iv,lt,level,shortname,longname,units,stdname)
                        
            for s in stationlist.keys():
                ncperdate[s].postproc_pz()
                ncperdate[s].postproc_tq()
                ncperdate[s].postproc_radiation()
                ncperdate[s].postproc_uv()
            print 'done with latest files'; #sys.exit(1)
