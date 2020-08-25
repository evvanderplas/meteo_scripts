#! /usr/bin/env python

import os,sys,time
import datetime as dt
import dateutil

def daterange(start_date, end_date = None, delta = (None,None)):
    '''
    makes an iterator of days within a period between startdate and enddate
    or uses a delta: eg delta = ('months',1)
    '''
    
    if end_date == None and delta == (None,None):
        print 'Provide either enddate or a timedelta (eg ("weeks",1))'
        sys.exit(1)
    #elif end_date == None and delta != (None,None):
    elif delta != (None,None):
        if delta[0] == 'hours':
            print 'yes hours yes'
            end_date = start_date + dateutil.relativedelta.relativedelta(hours=int(delta[1]))
            for n in range(int ((end_date - start_date).total_seconds()/3600.)):
                print start_date, dt.timedelta(hours=n)
                yield start_date + dt.timedelta(hours=n)
        elif delta[0] == 'days':
            end_date = start_date + dateutil.relativedelta.relativedelta(days=int(delta[1]))
        elif delta[0] == 'weeks':
            end_date = start_date + dateutil.relativedelta.relativedelta(weeks=int(delta[1]))
        elif delta[0] == 'months':
            end_date = start_date + dateutil.relativedelta.relativedelta(months=int(delta[1]))
        elif delta[0] == 'years':
            end_date = start_date + dateutil.relativedelta.relativedelta(years=int(delta[1]))

    for n in range(int ((end_date - start_date).days)):
        #print start_date, dt.timedelta(days=n)
        yield start_date + dt.timedelta(days=n)

def roundTime(dtm=None, roundTo=60):

   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """

   if dtm == None : dtm = dt.datetime.now()
   seconds = (dtm - dtm.min).seconds
   # // is a floor division, not a comment on following line:
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dtm + dt.timedelta(0,rounding-seconds,-dtm.microsecond)

def todtg(posdate):

    #print 'posix in ',posdate
    vdate      = dt.datetime.fromtimestamp(posdate)
    isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
    
    return vdate,isodate

def toposix(dtdate):

    #print 'datetime ',dtdate,
    if type(dtdate) == dt.datetime:
        posixdate  = int(time.mktime(dtdate.timetuple()))
    
    return posixdate

if __name__ == '__main__':

    print(dir(dateutil),dateutil.__file__)

    # test
    begindate = dt.datetime(2013,4,8,0)
    enddate   = dt.datetime(2013,4,10,0)
    
    # example of daterange:
    for single_date in daterange(begindate, enddate):
        print 'format, eg: ',single_date.strftime('%Y-%m_%d') #dt.datetime.strftime("%Y-%m-%d", single_date)
        for st in range(0,24,3):
            sttime =  single_date + dt.timedelta(hours = st)
            print 'with hours: ',sttime

    print 'test dateutil relativedelta '
    print begindate, begindate + dateutil.relativedelta.relativedelta(hours=1), begindate + dateutil.relativedelta.relativedelta(weeks=1), begindate + dateutil.relativedelta.relativedelta(months=1)
