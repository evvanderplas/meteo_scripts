#! /usr/bin/env python

import re,time,os,sys
import numpy as np
import datetime as dt
import netCDF4


def finddate(somefname):
    fc = re.compile('\S\D*(\d{4})(\d{2})(\d{2})(\d{2})\d*\D*(\d{3})\d*\D*$',re.VERBOSE)
    if 1:
        dtl = np.array(fc.search(somefname).groups(),dtype=int)
        #return dtl
        year,month,day,hour = dtl[0],dtl[1],dtl[2],dtl[3]
        leadtime = dtl[4]
        startdate = dt.datetime(year,month,day,hour,0,0) 
        validdate = dt.datetime(year,month,day,hour,0,0) + dt.timedelta(hours=int(leadtime))
        #print 'found ',startdate,validdate
        #writetime = int( date2num( validdate, vtime.units, calendar='standard') )
        return startdate,validdate

    else: #except:
        print 'Filename does not match yyyymmddhh bla lt bla',somefname
        return None,None

filenames = ['/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00000_GB', 
              '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00100_GB', 
              '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00200_GB', 
              '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00300_GB', 
              '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00400_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00500_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00600_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00700_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00800_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_00900_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01000_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01100_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01200_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01300_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01400_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01500_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01600_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01700_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01800_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_01900_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02000_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02100_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02200_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02300_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02400_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02500_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02600_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02700_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02800_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_02900_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03000_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03100_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03200_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03300_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03400_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03500_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03600_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03700_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03800_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_03900_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04000_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04100_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04200_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04300_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04400_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04500_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04600_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04700_GB', '/data/bens03/harmonie_data/GVDB/HARM_N25_201408120600_04800_GB']

tunits = 'seconds since 1970-01-01 00:00:00'
tstart = dt.datetime(1970,1,1,0,0)

for f in filenames:

    fpath,fname = os.path.split(f)

    refdate,validdate = finddate(fname)
    vtime_nc = netCDF4.date2num(validdate, tunits, calendar='standard') 
    vtime_mk = time.mktime(validdate.timetuple())
    vtime_dl = validdate - tstart

    #print 'netCDF4 date2num  : output, int(output)',validdate,vtime_nc, repr(vtime_nc),int(vtime_nc)
    #print 'time module mktime: output, int(output)',vtime_mk, int(vtime_mk)
    print 'time module mktime: output, int(output)',vtime_dl.total_seconds(),vtime_nc, vtime_mk, vtime_mk - vtime_nc, int(vtime_mk)
