#! /usr/bin/env python

import os,sys,time
import tarfile,sqlite3

import numpy as np
import datetime as dt

# convenience
from os.path import join as pjoin
from collections import OrderedDict as Ordict

def suf2(h):

    if int(h)<10: return '0'+str(h)
    elif 10 <= int(h) < 100: return str(h)
    elif 100<= int(h) < 1000: return '0'+str(int(h)/100)
    elif 1000<= int(h) < 10000: return str(int(h)/100)
    else:
        print 'wtf?', h; sys.exit(1)

def create_vfld_sqldb(dbfile,tablename = 'determ'):

    # connect to db
    conn   = sqlite3.connect(dbfile)
    c      = conn.cursor()

    #tablename  = tablename
    sql_create = 'create table {tab} (date integer, sdate integer, leadtime integer, stationid integer, mslp integer, precip real, clc real, t2m real, td2m real, relhum real, spechum real, tmax real, tmin real, wind real, winddir real, windmax real, gust real, gust6 real,  vis real, PRIMARY KEY (date,leadtime,stationid))'.format(tab = tablename)

    try:
        c.execute(sql_create)
        print 'Created table {tab} in database {db}'.format(tab = tablename, db = dbfile)
    except:
        #print 'Table {tab} already in database {db}'.format(tab = tablename, db = dbfile)
        pass

    conn.close()
    return 'Done'


def parse_act_ascii(vfile, verb = False,rexp = '^06+[23]'):

    import re
    stationpattern = re.compile(r'''{0} '''.format(rexp))

    if os.path.exists(vfile): # check if file exists

        # open file
        mo = open(vfile,'r')
        sqlvals = []
        stationdict = {}

        line = mo.readline()
        s = line.split(',')

        vdate = dt.datetime.strptime(s[1].strip()+s[2].strip(),'%Y%m%d%H%M') 
        #print vdate; sys.exit(1)

        for line in mo:
            s = line.split(',')
            #print s
            if re.match(rexp,s[0]) and s[5].strip() != '': 

                # test if there is already data from this station:
                if not stationdict.has_key(int(s[0])): stationdict[int(s[0])] = {}
                else:
                    #print 'inner dict already there'
                    pass

                if verb and s[4].strip() == 'TT': print 'T1.5m ',s
                

                if   s[4].strip() == 'ff': stationdict[int(s[0])]['wind'] = float(s[5]); print 'ff',s
                elif s[4].strip() == 'fx': stationdict[int(s[0])]['gust'] = float(s[5]); print 'fx',s
                elif s[4].strip() == 'tx': stationdict[int(s[0])]['t2m']= float(s[5]); print 't1.5x',s
                elif s[4].strip() == 'R1H': stationdict[int(s[0])]['precip']= float(s[5]); print 'pcp ',s

            #sys.exit(0)
        mo.close()
        if 0:
            sqlvals.append(
                #posixdate,isodate,h, # dat time group
                [station_id,          # station ID
                 mslp,precip,clc,     # pressure, precip, cloud cover
                 t2m,td2m,relhum, spechum, tmax, tmin, # temperatures etc
                 wind, winddir, windmax, gust, gust6, vis # wind etcetera
                 ]
                )

    for stat in  stationdict.keys(): 
        stationdict[stat]['date'] = vdate
        print stat, stationdict[stat]
    return stationdict



def parse_vfld_ascii(vfile, verb = False,rexp = '^06+[23]'):

    '''
    Routine that parses an ascii vfld file, uses the filename, returns values in right order

    can be shortened considerably with something like: (not (yet) done for readibility)

    indexlist = [0,9,4,10,11,7,13,14,15,8,6,5,18,16,17,12]
    sqlvals.append([ s[i] for in in indexlist])

    uses regular expression to discriminate in the stations from allsynop.list to be registered
    Here by default stations starting with 062?? and 063?? are taken (Netherlands)
    '''

    import re
    stationpattern = re.compile(r'''{0} '''.format(rexp))


    if os.path.exists(vfile): # check if file exists

        # open file
        mo = open(vfile,'r')
        sqlvals = []
        for line in mo:
            s = line.split()
            if len(s) == 19 and re.match(rexp,s[0]) : # hard-coded, not very elegant, but well, might be imported from, eh, fldextr.F ?
                #if verb: print s
                if verb and int(s[0]) == 6260: print 'De Bilt: ',s
                station_id = s[0]
                mslp = s[9]
                clc  = s[4]
                precip = s[10]
                spechum = s[11]
                t2m = s[7]
                td2m = s[13]
                tmax = s[14]
                tmin = s[15]
                relhum = s[8]
                wind = s[6]
                winddir = s[5]
                windmax = s[18]
                gust = s[16]
                gust6 = s[17]
                vis = s[12]

                # prepare data to be inserted with some sql statement analogous to creation syntax of db
                # sql = 'insert into {tab} values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'.format(tab = tablename)

                sqlvals.append(
                    #posixdate,isodate,h, # dat time group
                    [station_id,          # station ID
                    mslp,precip,clc,     # pressure, precip, cloud cover
                    t2m,td2m,relhum, spechum, tmax, tmin, # temperatures etc
                    wind, winddir, windmax, gust, gust6, vis # wind etcetera
                    ]
                    )

                if verb: print len(sqlvals), 'vals: ',sqlvals

                #cursor.execute(sql,sqlvals)

            else: # other info, not parsed
                if verb: print s, bool(re.match(rexp,s[0]))
                pass

        # clean up
        mo.close()

        return sqlvals

    else:
        print 'Not in untar: ',vfile
        return 0


def get_date(filename):

    '''

    '''
    import re
    if filename[0:4] == "vobs":
        ddt = re.compile('(\S*)(\d{10})(\D*)')
    elif filename[0:4] == "vfld":
        ddt = re.compile('(\S*)(\d{10})(\d{2})(\D*)')

    m = re.search(ddt,filename)
    if m:
        vdate = dt.datetime.strptime(m.group(2),'%Y%m%d%H')
        posixdate  = int(time.mktime(vdate.timetuple()))
        isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
        if filename[0:4] == 'vfld':
            h = int(m.group(3))
        else:
            h = 0
    else:
        raise ValueError('Could not find date from file ',vname)

    return posixdate,isodate,h

def parse_to_sqlite(vfiles,dbfile,tabname,vdate = None,leadtime = None,verbose = False):

    '''
    needs files as outputted by fldextr.f90 (vfiles), database and table name
    returns dbfile
    '''

    try:
        create_vfld_sqldb(dbfile,tablename = tabname)
    except:
        print 'Table exists: ',dbfile,tabname

    # connect to db
    conn   = sqlite3.connect(dbfile)
    c      = conn.cursor()


    if type(vfiles) != list:
        vfiles = [vfiles]

    # loop over vfld files:
    for vfile in vfiles:

        vpath,vname = os.path.split(vfile)
        vals        = parse_vfld_ascii(vfile, verb = False,rexp = '^06+[23]')
        if len(vals) == 0:
            if verbose: print 'no data in ',vfile
            continue #return None

        # decide which time to pick: the date given with the call of the function, or the date in the filename
        if vdate is not None:
            posixdate  = int(time.mktime(vdate.timetuple()))
            isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
        else:
            posixdate,isodate,h = get_date(vname)

        if vname[0:4] == 'vfld':   # check if combination of date, leadtime, station are already there
            h = int(vname[-2:]); print h, isodate; #return None#sys.exit(1)
            c.execute("SELECT count(*) FROM {tab} WHERE date = ? AND leadtime = ? AND  stationid = ? ".format(tab=tabname), (posixdate,h,vals[0][0]))
        elif vname[0:4] == 'vobs': # check if combination of date, station are already there
            c.execute("SELECT count(*) FROM {tab} WHERE date = ? AND  stationid = ? ".format(tab=tabname), (posixdate,vals[0][0]))

        data=c.fetchone()[0];print 'count of data: ',data
        if data==0:
            print 'new, insert:',
            if verbose: print('There is no date {}, leadtime {} for station {}'.format(posixdate,h,vals[0][3]))
            pass
        else:
            print 'already there, moving on'
            if verbose: print('Component %s found in %s row(s)'%(posixdate,data))
            continue #break


        # times not necessarily available within file
        for row in vals:

            row.insert(0,posixdate)
            row.insert(1,isodate)
            row.insert(2,h)

            sql = 'insert into {tab} values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'.format(tab = tabname)
            #if verbose: print sql,row
            try:
                c.execute(sql,row)
            except: # sqlite3.IntegrityError:
                print('Primary key already in table ',row)
                continue

    # end of loop over files

    print 'committing to database',dbfile
    conn.commit()

    # close connection to database
    conn.close()

    return dbfile


def parse_dir(tardir,case = None,casetar= None,tempdir = './temp',
              latest = True, date = None, hbefore = 24,hafter=24, hcycle = 3,
              verbose = False):

    '''
    Uses a cursor of the (sqlite?) database connection
    (not having to open and close the database file for every single parsed file)
    e.g. :
    connection = sqlite3.connect(dbfile)
    cursor     = connection.cursor()
    '''

    maxlead = 48 # max leadtime

    parsedir = tempdir #'./temp'
    try:
        os.mkdir(parsedir)
    except:
        print 'Directory already in place: ', parsedir


    ## If latest = True only the files from the latest run are parsed
    ## else the files of the date until hbefore hours before the date are parsed

    #import checkdate
    if latest:
        date = dt.datetime.today()
        hours = range(-hbefore,max(hafter,hcycle),hcycle); print 'For ',dtg,hours

    elif date is not None:
        if type(date) in (int,str):
            try:
                date = dt.datetime.strptime(str(date),'%Y%m%d%H')
            except:
                print "could not derive starting time from ",date
                return None

        # Parse all available files for a previous date
        #hours = range(0,int(dtg[8:10]),3); print 'For ',dtg,hours
        hours = range(-hbefore,hafter,hcycle); print 'For ',dtg,hours

    else:
        print 'No date given: try latest = True for latest run, or a (recent) date'
        return None

    # just to make sure
    date.replace(hour=date.hour - date.hour%hcycle,minute=0,second=0,microsecond=0 )
    dtg  = date.strftime('%Y%m%d%H')

    ## The name of the tarfile might be different from the vfld-files it contains
    ## case is for the name of the tarfile
    ## casetar for the vfld-files contained in the tarred directory
    if casetar == None: casetar = case


    dbdir   = './'
    dbfile  = os.path.join(dbdir,'vfld_'+casetar+'_'+dtg[0:6]+'.db')

    # call create
    tablename = 'determ'

    # is file/table already there?
    try:
        create_vfld_sqldb(dbfile,tablename = tablename)
    except:
        print 'Table exists'

    # connect to db
    conn   = sqlite3.connect(dbfile)
    c      = conn.cursor()

    ## create vfld tarfile based on dtg and case
    for st in hours: # try to take latest/newest version of this file
        dtg2 = (date + dt.timedelta(hours = st)).strftime('%Y%m%d%H')
        tfile  = os.path.join(tardir,'vfld'+case+dtg2[0:10]+'.tar.gz')
        if not os.path.exists(tfile):
            print 'Does not exist: ',tfile
            pass
        else:
            print st, tfile, os.path.exists(tfile)

            untar_vfld(tfile,tempdir = parsedir,verb = True)

            # take starting date-time combined with leadtime as key
            #stdate     = dt.datetime.strptime(dtg[0:8]+suf2(st),'%Y%m%d%H')  #print stdate; sys.exit(1)
            stdate     = dt.datetime.strptime(dtg2[0:10],'%Y%m%d%H')  #print stdate; sys.exit(1)
            posixdate  = int(time.mktime(stdate.timetuple()))
            isodate    = int(dtg[0:10]+'0000')

            for h in range(maxlead+1): # h is leadtime

                ## create vfld filename based on dtg and casetar
                vfile = os.path.join(parsedir,'vfld'+casetar+dtg2[0:10]+suf2(h))
                if not os.path.exists(vfile):
                    print 'Not in tarfile',h,vfile
                    continue #break

                ## parse vfld file and check if record already exists in database:
                vfldvals = parse_vfld_ascii(vfile,verb=verbose) # rexp = '^06+[23]' for dutch stations
                print 'station id to check:',vfldvals[0][0]

                c.execute("SELECT count(*) FROM determ WHERE date = ? AND leadtime = ? AND  stationid = ? ", (posixdate,h,vfldvals[0][0]))
                data=c.fetchone()[0];print 'count of data: ',data,
                if data==0:
                    print 'new, insert'
                    if verbose: print('There is no date {}, leadtime {} for station {}'.format(posixdate,h,vfldvals[0][3]))
                    pass
                else:
                    print 'already there, moving on'
                    if verbose: print('Component %s found in %s row(s)'%(posixdate,data))
                    continue #break


                ## make an array/list of data to be inserted, and pass it to sql statement:
                #  times not necessarily available within file
                for row in vfldvals:

                    row.insert(0,posixdate)
                    row.insert(1,isodate)
                    row.insert(2,h)

                    sql = 'insert into {tab} values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'.format(tab = tablename)
                    c.execute(sql,row)

            if latest: break # take only last one

    print 'committing to database',dbfile
    conn.commit()

    # close connection to database
    conn.close()

    return dbfile


def allsynop_to_sqlite(synopfile,dbfile,tabname= 'stations',verbose = False):
    
    '''
    Put the file with all the description data per synop station (allsynop.list)
    as a separate table in a database
    '''

    # connect to db (use sqlite3)
    conn   = sqlite3.connect(dbfile)
    conn.text_factory = str
    c      = conn.cursor()

    sql    = "create table {tab} (sid integer, lat real, lon real, elev real, name text)".format(tab = tabname)
    try:
        if verbose: print sql
        c.execute(sql)
    except sqlite3.OperationalError:
        if verbose: print 'Table exists: ',tabname
        return dbfile
        #pass

    # parse synop stations
    sf = open(synopfile, 'r')
    sf.readline()
    for line in sf:
        try:
            station = line.split()
            sid = station[0]
            lat = station[1]
            lon = station[2]
            elev = station[3]
            name = ' '.join(l for l in station[4:])

            sql    = "insert into  {tab} values (?,?,?,?,?)".format(tab = tabname)
            if verbose: print sql, sid,lat,lon,elev,name
            c.execute(sql,(sid,lat,lon,elev,name))
        except:
            pass

    print 'committing to database',dbfile,', closing'
    conn.commit()

    # close connection to database
    conn.close()

    return dbfile

##---------------------------------------------------------------
    '''
                station_id = s[0]
                mslp = s[9]
                clc  = s[4]
                precip = s[10]
                spechum = s[11]
                t2m = s[7]
                td2m = s[13]
                tmax = s[14]
                tmin = s[15]
                relhum = s[8]
                wind = s[6]
                winddir = s[5]
                windmax = s[18]
                gust = s[16]
                gust6 = s[17]
                vis = s[12]
     '''

parameters = Ordict({
    'startdate':{
            'shortname':'startdate',},
    'validdate':{
            'shortname':'validdate',},
    'leadtime':{
            'shortname':'leadtime',},
    'sid'  :{
            'shortname':'sid',},
    'lat'  :{
            'shortname':'lat',},
    'lon'  :{
            'shortname':'lon',},
    'z'  :{
            'shortname':'z',},
    # to be able to have the header info in the same table
    'NN'  : {'name':'Cloud Cover',
             'shortname':'CC',
             'pos':4},
    'DD'  : {'name':'Wind direction 10m',
             'shortname':'D10m',
             'pos':5},
    'FF'  : {'name':'Wind speed 10m',
             'shortname':'FF10m',
             'pos':6},
    'TT'  : {'name':'Temperature 2m',
             'shortname':'T2m',
             'pos':7},
    'TD' : {'name':'Dew point temperature 2m',
            'shortname':'TD2m',
            'pos':13},
    'RH'  : {'name':'Relative humidity 2m',
             'shortname':'RH2m',
             'pos':8},
    'QQ'  : {'name':'Specific humidity 2m',
             'shortname':'Q2m',
             'pos':11},
    'PSS' : {'name':'Surface pressure',
             'shortname':'PMSL',
             'pos':9},
    'PS'  : {'name':'Surface pressure',
            'shortname':'SP',
            'pos':9},
    'VI' : {'name':'Visibility',
            'shortname':'vis',
            'pos':12},
    'TX' : {'name':'Maximum temperature 2m',
            'shortname':'Tmax',
            'pos':14},
    'TM' : {'name':'Minimum temperature 2m',
            'shortname':'Tmin',
            'pos':15},
    'GW' : {'name':'Wind gust 10m',
            'shortname':'FX10m',
            'pos':16},
    'GM' : {'name':'Max wind gust 10m 6 h',
            'shortname':'FX10m6H',
            'pos':17},
    'WX' : {'name':'Max wind 10m 6 h',
            'shortname':'FF10m6H',
            'pos':17},
    'PE1' : {'name':'Precipitation, accum 1 h',
            'shortname':'pcpacc1',
            'pos':17},
    'PE3' : {'name':'Precipitation, accum 3 h',
            'shortname':'pcpacc3',
            'pos':17},
    'PE6' : {'name':'Precipitation, accum 6 h',
            'shortname':'pcpacc6',
            'pos':17},
    'PE' : {'name':'Precipitation, accum 12 h',
            'shortname':'pcpacc12',
            'pos':17},
    'PE24' : {'name':'Precipitation, accum 24 h',
            'shortname':'pcpacc24',
            'pos':17},
    'PP'  : {'name':'pressure (PL)',
             'shortname':'P',
             'pos':9},
    'FI'  : {'name':'Geopotential height',
            'shortname':'GP',
            'pos':9},
   })
vobsmaxlist   = ['NN','FF','DD','TT','TD','RH','QQ','PSS','PS','VI','TX','TM','GW','GM','WX','PE1','PE3','PE6','PE','PE24']
vobsmaxlistv2 = ['NN','DD','FF','TT','RH','PS','PE','QQ','VI','TD','TX','TM','GW','GM','WX']
'''
v2:
synop    write(16,999)stnr_s(j),lat_s(j),lon_s(j),fi_s(j),
     x                nn_s(j),dd_s(j),ff_s(j),tt_s(j),
     x                rh_s(j),pp_s(j),rr_s(j),qq_s(j),
     x                vis_s(j),td_s(j),
     $                TX_s(j),TM_s(j),GU_s(j),G6_s(j),WX_s(j)

temp        write(16,1001)pres_t(jlev),fi_t(j,jlev),
     x                    tt_t(j,jlev),rh_t(j,jlev),
     x                    dd_t(j,jlev),ff_t(j,jlev),
     x                    qq_t(j,jlev),td_t(j,jlev)

'''
tempmaxlist   = ['PP', 'FI', 'TT', 'RH', 'DD', 'FF', 'QQ', 'TD']
tempmaxlistv2 = ['PP', 'FI', 'TT', 'RH', 'DD', 'FF', 'QQ', 'TD']

def toposix(dtdate):

    #print 'datetime ',dtdate,
    if type(dtdate) == dt.datetime:
        posixdate  = int(time.mktime(dtdate.timetuple()))
    
    return posixdate

def todtg(posdate):

    #print 'posix in ',posdate
    vdate      = dt.datetime.fromtimestamp(posdate)
    isodate    = int(vdate.strftime('%Y%m%d%H%M')+'00')
    
    return vdate,isodate

def untar_vfld(tfile,tempdir = './temp',verb = False):

    import tarfile

    f = tarfile.open(tfile,'r:*')

    tarlist = []
    fm = f.getmembers()
    for m in fm:
        if verb:
            print m.name
        tarlist.append(m.name)

    f.extractall(path = tempdir)

    f.close()
    return tarlist #'Done'


class ascii:

    def __init__(self,asciifile,date=None):
        self.file = asciifile
        self.dtg  = asciifile[4:]
        if date == None:

            try:
                self.date = dt.datetime.strptime(self.dtg,'%Y%m%d%H')
                print self.date
            except:
                print 'Hm, date?',asciifile

def process_v(vfile,outdir = './',tempdir= None,model=None,vdate=None):

    print os.path.splitext(vfile)
    fpath,fname = os.path.split(vfile)
    if fname[0:4] == 'vobs':
        leadtime = 0
        
    if os.path.splitext(vfile)[1] == '.gz':
    
        if tempdir == None:
            tempdir= pjoin(outdir,'temp')
            try:
                os.makedirs(tempdir)
            except:
                print 'temp in place'

        tempdir = os.path.abspath(tempdir)
        print 'Extracting to ',tempdir

        tlist = untar_vfld(vfile,tempdir = tempdir)
        for t in tlist:
            tpath,tname = os.path.split(t)
            tfile = pjoin(tempdir,tpath,tname)
            if tname[0:4] == 'vobs':
                if vdate == None:
                    vdate    = dt.datetime.strptime(tname[4:],'%Y%m%d%H')
            elif tname[0:4] == 'vfld':
                leadtime = int(tname[-2:])
                if model == None:
                    model = tname[4:-12]
                    print 'Distilled modelname:',model
                if vdate == None:
                    vdate = dt.datetime.strptime(tname[-12:-2],'%Y%m%d%H')
            else:
                print 'Not sure how to parse:',t
                sys.exit(0)

            parse_v(tfile,vdate=vdate,model=model,leadtime=leadtime,outdir = outdir)

            print 'Removing ',tfile
            os.remove(tfile)
    
    else:
        if fname[0:4] == 'vobs':
            leadtime = 0
            vdate    = dt.datetime.strptime(fname[4:],'%Y%m%d%H')
        elif fname[0:4] == 'vfld':
            modelname = fname[4:-12]
            vdate     = dt.datetime.strptime(fname[-12:2],'%Y%m%d%H')
            leadtime  = int(fname[-2:])
        
        parse_v(vfile,vdate=vdate,model=model,leadtime=leadtime,outdir = outdir)


def parse_v(vobsfile,vdate=None,model='default',leadtime=0,outdir = './'):

    ## Some boilerplate code, sanity check:

    ## Check based on file name: ok for now? 
    fpath,fname = os.path.split(vobsfile)
    if   fname[0:4] == 'vobs': vtype = 'obs'
    elif fname[0:4] == 'vfld': vtype = 'fc'

    if not os.path.exists(vobsfile):
        print 'File does not exist',vobsfile
        return -1

    try:
        f = open(vobsfile,'r')
        fl = f.readline()
        print fl
        #f.close()
    except:
        print 'File not legible:',vobsfile
        return -1,-1

    try:
        os.makedirs(outdir)
    except:
        print 'Outdir in place',outdir

    if vdate == None:
        try:
            vdate = dt.datetime.strptime(fname[4:],'%Y%m%d%H')
            print vdate
        except:
            print 'Hm, date?',vobsfile
            sys.exit(0)

    svals,tvals = [],[]
    if int(fl.split()[2]) in (2,3):
        print 'Processing type 2',vdate,leadtime
        spars,svals,tpars,tvals = parse_v2(f,toposix(vdate),leadtime)
        #print 'after v2',spars
        #print 'after v2',tpars
    elif fl.split()[2] == '4':
        print 'Processing type 4',vdate,leadtime
        spars,svals,tpars,tvals = parse_v4(f,toposix(vdate),leadtime) 
        #print 'after v4',spars
        #print 'after v4',tpars

    if svals != []:
        if  vtype == 'obs':
            dbfile = 'synop_{ym}.db'.format(ym = vdate.strftime('%Y%m'))
        elif vtype == 'fc':
            dbfile = 'synop_{m}_{ym}.db'.format(m=model,ym = vdate.strftime('%Y%m'))
        create_vobs_db(dbfile,outdir = outdir)
        print "inserting synop"
        insert_vobs_db(dbfile,svals,data='S',parlist = spars,outdir = outdir)
        print "inserting temp",len(tvals),tpars
        insert_vobs_db(dbfile,tvals,data='T',parlist = tpars,outdir = outdir)
    else:
        print 'No parseable data in ',vobsfile

    # cleanup
    f.close()

    return svals,tvals #'done'

def parse_v2(fobj,startdate,leadtime):

    '''
    Parsing a type 2 vobs file
    Return a list of lists to insert into SQLite table
    '''

    # create a human readable date as validdate:
    vdate,validdate = todtg(startdate + leadtime)

    synopheader = ['sid','lat','lon','z']
    vmlist = ['startdate','validdate','leadtime'] + synopheader + vobsmaxlist
    v2list = ['startdate','validdate','leadtime'] + synopheader + vobsmaxlistv2

    tmlist = ['startdate','validdate','leadtime'] + ['sid'] + tempmaxlist
    t2list = ['startdate','validdate','leadtime'] + ['sid'] + tempmaxlistv2

    #reorder  = [vmlist.index(n) for n in v2list]
    #invorder = [v2list.index(n) for n in vmlist]
    #vobsar    = np.array(vmlist)
    #vobsarinv = np.array(v2list)

    spars,tpars = v2list,t2list
    svals,tvals = [],[]

    for l in fobj:
        s = l.split()
        if len(s) == 19:
            a = np.array([startdate,validdate,leadtime]+s,dtype=float)
            svals.append(a)
        elif len(s) == 4:
            sid = s[0]
            #print sid
        elif len(s) == 8: #9:
            a = np.array([startdate,validdate,leadtime]+[sid]+s,dtype=float)
            tvals.append(a)

    return spars,svals,tpars,tvals

def parse_v4(fobj,startdate,leadtime):

    '''
    Parsing a type 4 vobs file
    Return a list of lists to insert into SQLite table
    '''

    # create a human readable date as validdate:
    vdate,validdate = todtg(startdate + leadtime)

    synopheader = ['startdate','validdate','leadtime','sid','lat','lon','z']
    tempheader  = ['startdate','validdate','leadtime','sid']

    ## NB EvdP hard code:
    #startdate   = 139009882772
    
    fobj.seek(0)
    fl = fobj.readline()
    sl = fobj.readline()
    nrofsynop = int(fl.split()[0])
    nrofspars = int(sl.split()[0])
    nroftpars = -1
    #print fl,sl; return 0

    fobj.seek(0)
    readsynop,readpl = False, False

    spars,tpars = [],[]
    svals,tvals = [],[]
    
    for l in fobj:
        s = l.split()
        #print 'v4',startdate,validdate,leadtime,len(spars),len(synopheader+spars),len(s),nrofspars
        #print 'synop length',len(s),len(synopheader)-3 + nrofspars

        # do the parsing logic...
        if len(s) == 1 and readsynop:
            if not readpl:
                nrofpl = int(s[0]) 
                readpl = True
            else: # if nrofpl has been read
                nroftpars = int(s[0])
        elif len(s) == 2 and s[0][0].isalpha():
            if not readsynop:
                spars.append(s[0])
                #print spars
            else:
                tpars.append(s[0])
                #print tpars
        elif len(s) == len(synopheader)-3 + nrofspars: #len(synopheader + pars):
            #print 'reading synop data',s[0]
            readsynop = True # once a synop line is found, reset: afterwards come temps
            svals.append([startdate,validdate,leadtime]+s)
        elif len(s) == len(synopheader)-3:
            sid = s[0]
        elif len(s) == nroftpars: #len(synopheader + pars):
            #print 'reading temp data',s[0]
            tvals.append([startdate,validdate,leadtime]+[sid]+s)
    
    return synopheader+spars,svals,tempheader+tpars,tvals

def sqlexec(cursor,sql):

    # execute simple with try
    try:
        cursor.execute(sql)
        return 1
    except sqlite3.OperationalError:
        print 'Nope: ',sql
        return 0

def sqlinsert(cursor,sql,valslist):

    # execute many with try
    if 1: #try:
        cursor.executemany(sql,valslist)
        return 1
    else: #except sqlite3.OperationalError:
        print 'Nope: ',sql
        return 0

def create_vobs_db(dbname,outdir = './'):

    dbfile = pjoin(outdir,dbname)
    if os.path.exists(dbfile): return 0

    sqlcreate_hs = 'CREATE TABLE synop ('
    sqlcreate_ht = 'CREATE TABLE temp ('
    header_s     = 'startdate integer, validdate integer, leadtime integer, sid integer, lat real, lon real, z real, '
    header_t     = 'startdate integer, validdate integer, leadtime integer, sid integer, '
    body_s       = ', '.join([parameters[k]['shortname'] + ' real' for k in vobsmaxlist])
    body_t       = ', '.join([parameters[k]['shortname'] + ' real' for k in tempmaxlist])
    sqlcreate_es = ', PRIMARY KEY (sid,leadtime,startdate))'
    sqlcreate_et = ', PRIMARY KEY (sid,leadtime,startdate,P))'
    sqls = sqlcreate_hs + header_s + body_s + sqlcreate_es
    sqlt = sqlcreate_ht + header_t + body_t + sqlcreate_et
    print sqls; #return 0

    # connect to db
    conn   = sqlite3.connect(dbfile)
    c      = conn.cursor()

    s = sqlexec(c,sqls)
    t = sqlexec(c,sqlt)

    conn.close()
    return 'Done'

def insert_vobs_db(dbname,valuelist,parlist = None, data='S',outdir='./'):

    dbfile = pjoin(outdir,dbname)
    if not os.path.exists(dbfile): 
        print 'DB does not exist (yet):', dbfile
        return 0
    if len(valuelist) > 1:
        nr = len(valuelist[0])
        print 'Nr of values:',nr
    else:
        print 'values, eh',len(valuelist),valuelist
        return 'Done (nothing)'
        #sys.exit(0)

    tables= {'S':'synop','synop':'synop',
             'T':'temp','temp':'temp',
             }
    
    if parlist == None:
        valph    = '('+','.join('?' for i in range(nr)) + ')'
        sqlbulk  = 'INSERT into {tab} VALUES {ph}'.format(tab=tables[data],ph=valph)
    else:  #parlist != None:
        #print parlist
        inslist  = '('+','.join([parameters[p]['shortname'] for p in parlist]) + ')'
        #print inslist
        valph    = '('+','.join('?' for i in range(nr)) + ')'
        sqlbulk  = 'INSERT into {tab} {cols} VALUES {ph}'.format(tab=tables[data],cols=inslist,ph=valph)

    #print valuelist[0]
    #print valuelist[-1]
    #print sqlbulk
    lts,sids = [],[]
    for l in valuelist:
        lts.append(l[2])
        sids.append(l[3])

    #print sids
    print 'length synop all :',len(sids)
    sidset = set(sids)
    ltset = set(lts)
    print 'length synop unique',len(sidset),len(ltset)
    #sys.exit(0)

    # connect to db
    conn   = sqlite3.connect(dbfile)
    c      = conn.cursor()

    try:
        ins = sqlinsert(c,sqlbulk,valuelist)
        conn.commit()
    except sqlite3.IntegrityError:
        print 'Columns not unique',valuelist[0]
        pass

    conn.close()
    return 'Done'
    

if __name__ == '__main__':


    # For (unit) testing purposes!
    print '********* testing the parsing of vfld files to sqlite database **********'


    today = dt.datetime.today()
    today = dt.datetime(2014,3,25,18)

    dtgh  = dt.datetime.strftime(today,'%Y%m%d%H')
    testdir = '/usr/people/plas/python/tools/vobs'
    #testdir = '/home/plas/python/tools/vobs'
    vobsfile = 'vobs{dtgh}'.format(dtgh=dtgh)
    
    vobs2 = pjoin(testdir,'v2',vobsfile)
    vobs4 = pjoin(testdir,'v4',vobsfile)
    print vobs2, os.path.exists(vobs2)
    print vobs4, os.path.exists(vobs4)

    if 0:
        s,t = parse_v(vobs2)
        #1395766800  |6240|52.32|4.79|-3.0|0.0|6.0|30.0|280.8|274.1|62.5   |0.004016 |      |1012.2|30000.0|284.0|-99.0|8.0|11.0|7.0|     |     |     |0.0|
        synopvals,tempvals = parse_v(vobs4)
        #139009882772|6240|52.32|4.79|-3.0|0.0|6.0|30.0|280.8|274.1|62.4686|0.0040159|1012.7|1012.2|30000.0|284.0|-99.0|8.0|11.0|7.0|-99.0|-99.0|-99.0|0.0|-99.0    

    if 0:
        today = dt.datetime(2014,4,5,18)
        model = 'racmoturbzch02'
        for h in range(25):
            vfldfile = 'vfld{m}{ymdh}{lt}'.format(m=model,ymdh=today.strftime('%Y%m%d%H'),lt=str(h).zfill(2))
            vfld  = pjoin(testdir,'vfld',vfldfile)
            synopvals,tempvals = parse_v(vfld,vdate=today,leadtime=h*3600)

    if 1:
        today = dt.datetime(2014,4,7,0)
        model = 'racmoturbzch02'
        vtarfile  = 'vfld{m}{ymdh}.tar.gz'.format(m=model,ymdh=today.strftime('%Y%m%d%H'))
        vfld = pjoin(testdir,'vfld',vtarfile)
        process_v(vfld)
        #tfiles = untar_vfld(vfld,tempdir = pjoin(testdir,'vfld','temp'),verb = True)

    sys.exit(0)

    parsedir = './temp'
    try:
        os.mkdir(parsedir)
    except:
        print 'Directory already in place: ', parsedir

    tardir = '/net/bens03/harmonie_data/GVDB'
    tardir = '/net/bhw284/nobackup/users/tijm/verg/HARM36'

    case    = '36h14'
    maxlead = 48


    vfld_dbfile = parse_dir(tardir,case = case, tempdir = './temp',latest = True, 
                            date = today - dt.timedelta(hours = 48),
                            hbefore = 9,hafter = 0, hcycle=3,
                            verbose = False)

    print 'Parsed ',tardir,'\nto ',vfld_dbfile
