#! /usr/bin/env python

import os,sys
import datetime as dt
import netCDF4 as nc

def get_obs(date1,station_id = 260, date2 = None,df=False,verb = False):
    
    '''
    Get observations from KIS database using the cgi interface
    input: startdate and enddate (date1, date2 is optional, defaults to one day)
    output in lists (default: date,station id, t2m, precipitation, wind, pressure) or as 
    a pandas DataFrame (if df == 1/True)

    '''


    import urllib, urllib2

    if type(date1) == dt.datetime:

        datum1 = dt.datetime.strftime(date1,'%Y%m%d%H%M')

        if date2 == None:
            date2  = date1+dt.timedelta(hours=23)

        datum2 = dt.datetime.strftime(date2,'%Y%m%d%H%M');print 'until ',datum2

        vars="T:RH:FF:FH:FX:U:N"

        form_data = \
        {"byear" :date1.strftime('%Y'), #datum1[0:4],
         "eyear" :date2.strftime('%Y'), #datum2[0:4],
         "bmonth":date1.strftime('%m'), #datum1[4:6],
         "emonth":date2.strftime('%m'), #datum2[4:6],
         "bday"  :date1.strftime('%d'), #datum1[6:8],
         "eday"  :date2.strftime('%d'), #datum2[6:8],
         "bhour":'0',  #datum1[8:10],
         "ehour":'23', #datum2[8:10],
         'vars'  : vars, # optional
         #"stations":str(station_id), 
         "lang":"nl",
         }

        ## here check if station is integer or a list:
        if type(station_id) == int:
            form_data["stations"] = str(station_id)
        elif type(station_id) == list:
            form_data = form_data.items() # no longer a dict!
            stations = [('stations',str(s)) for s in station_id]
            form_data.extend(stations)
            
        print form_data
        ## mold it into url and request from website:
        url = 'http://www.knmi.nl/klimatologie/uurgegevens/getdata_uur.cgi'
        params = urllib.urlencode(form_data)
        if verb: print 'url,params = ',url,params

        response = urllib2.urlopen(url, params)
        lines    = response.readlines()

        ## parse/process output
        date,sid,t2m,pcp,ff,fh,fx,clc = [],[],[],[],[],[],[],[]
        for line in lines:
            if line[0] != '#':
                # this is a data line
                #print len(line),line
                #values = [float(v) for v in line.split(',') if v.strip().isdigit()]
                # ordering is important:...

                #values = [0 if not v.strip().isdigit() else float(v) for v in line.split(',')]
                try:
                    values = [float(v) for v in line.split(',')]
                except:
                    print 'Strange value in line:',line
                    print 'next'
                    continue

                dtg = dt.datetime.strptime(str(int(values[1])),'%Y%m%d') + dt.timedelta(hours = int(values[2]))
                date.append(dtg)
                sid.append(int(values[0]))
                t2m.append(values[3]/10.)
                if values[4]>=0:
                    pcp.append(values[4]/10)
                else:
                    pcp.append(0.)
                ff.append(values[5]/10)
                fh.append(values[6]/10)
                fx.append(values[7]/10)
                if values[9]>=0: #N
                    clc.append(values[9])
                else:
                    clc.append(0.)

    #for it in range(len(date)): print date[it],t2m[it]
    #sys.exit(0)

    if df:

        dates = [int(d.strftime('%Y%m%d%H')) for d in date] # make datetime objects
        sids  = [int(s+6000) for s in sid] # international station ID
        t2m   = [t + 273.15 for t in t2m]  # Kelvin au lieu de Celsius
        wind  = [f for f in ff]  # wind
        apcp  = [p for p in pcp]  # precip

        try:
            import pandas as pd
        except:
            print 'pandas not available, returning lists instead'
            return dates,sids,t2m,pcp,ff,fh,fx,clc

        obsdf = pd.DataFrame({'validdate':dates,
                                'sid':sids,
                                't2m':eval('t2m'),
                                'wind':eval('wind'),
                                'apcp':eval('apcp')
                            })
        return obsdf

    return date,sid,t2m,pcp,ff,fh,fx,clc

'''
HTTP:
http://birdexp03.knmi.nl/plieger/data/kmds/10M/20131217/

WMS:
http://birdexp03/cgi-bin/plieger/kmds.cgi? 
'''

def hour_range(start_date, end_date):
    for n in range(int ((end_date - start_date).total_seconds()/3600 ) ):
        yield start_date + dt.timedelta(hours=n)

def read10min(ncfile):

    uurdat = nc.Dataset(mync,'r')
    times  = [dt.datetime.fromtimestamp(posdate) for posdate in uurdat.variables['time'][:]]
    print times
    stations  = uurdat.variables['station'][:]
    t2m   = uurdat.variables['ta'][:] # temperature
    td2m  = uurdat.variables['td'][:] # dew point t
    ff    = uurdat.variables['ff'][:] # wind speed
    fx    = uurdat.variables['fx'][:] # wind gust
    rg    = uurdat.variables['rg'][:] # precip intens (mm/h)
    rh    = uurdat.variables['rh'][:] # relative humidity
    nc    = uurdat.variables['nc'][:] # total cloud cover
    print fx[stations == 'A260a']


def parse_10min(begindate, enddate = None,df = False):

    '''
    read 10 min files from KMDS site
    '''
    
    import netCDF4 as nc
    import numpy as np
    import urllib, urllib2
    import daterange

    if type(begindate) == dt.datetime: 
        begindate = begindate.replace(hour=0, minute=0, second=0, microsecond=0)
    if enddate == None: enddate = dt.datetime.today()

    #statids = 

    print begindate, enddate
    for t in hour_range(begindate,enddate): #daterange.daterange(begindate,enddate,delta=('hours',1)):
        print t, t.hour, str(t.hour).zfill(2)

        beginh,endh = t, t + dt.timedelta(seconds=3600)
        url = 'http://birdexp03.knmi.nl/plieger/data/kmds/10M/{dtgs}/KMDS__TEST_P___10M_OBS_L2__{dtgs}T{hb}0000_{dtge}T{he}0000_0001.nc'.format(dtgs = t.strftime('%Y%m%d'),hb=str(t.hour).zfill(2),dtge=endh.strftime('%Y%m%d'),he=str(endh.hour).zfill(2))
        print url
        mync = 'test_{h}.nc'.format(h=t.strftime('%Y%m%d%H'))

        if not os.path.exists(mync):
            print 'Downloading file ',mync
            urllib.urlretrieve(url,mync)
            print mync, os.path.exists(mync)
        # response = urllib2.urlopen(url)
        try:
            uurdat = nc.Dataset(mync,'r')
        except RuntimeError:
            print 'Trouble opening ',mync
            continue

        times  = np.array([dt.datetime.fromtimestamp(posdate) for posdate in uurdat.variables['time'][:]])
        print times
        stations  = uurdat.variables['station'][:]
        t2m   = uurdat.variables['ta'][:] # temperature
        td2m  = uurdat.variables['td'][:] # dew point t
        ff    = uurdat.variables['ff'][:] # wind speed
        fx    = uurdat.variables['fx'][:] # wind gust
        rg    = uurdat.variables['rg'][:] # precip intens (mm/h)
        rh    = uurdat.variables['rh'][:] # relative humidity
        tcc    = uurdat.variables['nc'][:] # total cloud cover
        print fx[stations == 'A260a']
        uurdat.close()

        try:
            timesday = np.concatenate([timesday,times])
            t2mday   = np.concatenate((t2mday,t2m),axis =1)
            td2mday  = np.concatenate((td2mday,td2m),axis =1)
            ffday    = np.concatenate((ffday,ff),axis =1)
            fxday    = np.concatenate((fxday,fx),axis =1)
            rgday    = np.concatenate((rgday,rg),axis =1)
            rhday    = np.concatenate((rhday,rh),axis =1)
            ncday    = np.concatenate((ncday,tcc),axis =1)
        except NameError:
            timesday,t2mday,td2mday,ffday,fxday,rgday,rhday,ncday = times,t2m,td2m,ff,fx,rg,rh,tcc

    if df: # not yet
        try:
            import pandas as pd
            return pd.Dataframe()
        except ImportError:
            print 'Pandas not available'
            return None
    else:
        return {'times':timesday, 
                'stations':stations, 
                't2m':t2mday, 
                'td2m':td2mday, 
                'ff':ffday, 
                'fx':fxday, 
                'apcp':rgday,
                'rh':rhday,
                'tcc':ncday 
                }



if __name__ == '__main__':

    
    tenmindict = parse_10min(dt.datetime.today(), enddate = None,df=False)
    sys.exit(0)

    adate = dt.datetime(2013,10,13)
    sodate,sid,t2m,pcp,ff,fh,fx,clc = get_obs(adate,verb=True)
    print pcp
    sodate,sid,t2m,pcp,ff,fh,fx,clc = get_obs(adate,station_id = [240,260],verb=True)
    

