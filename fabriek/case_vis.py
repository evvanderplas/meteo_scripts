#! /usr/bin/env python

import os,sys,glob,time,numpy
import datetime as dt
import sqlite3

backend = 'Agg'
import matplotlib
matplotlib.use(backend)

import matplotlib.pyplot,pylab
import matplotlib.cm as CM

import sources
import tools.sql_ts_vis as vis
import tools.select as select

from verif.tools import parse_vfld

def mkdir_p(path):
    try:
        os.makedirs(path)
    except:
        print 'Outdir already there:',path
    return path


if __name__ == '__main__':

    
    # some constants, list of stations:
    sid = (6260,6235,6370,6380,6240,6251,6270,6280)
    #sid = (6260,)
    hh = 36
    hb = 24
    hc = 3

    topdir  = '/net/bhw379/nobackup_1/users/plas/fabriek'
    tempdir = '/nobackup_1/users/plas/temp'
    filelist = {}

    ydate = dt.datetime.today() - dt.timedelta(hours = 24)
    date = dt.datetime(ydate.year, ydate.month, ydate.day); #print date; sys.exit(1)
    date = dt.datetime(2013,10,28)
    date = dt.datetime(2013,11,11)
    date = dt.datetime(2013,12,5)
    dtg  = date.strftime('%Y%m%d%H')


    #for fcsource in (sources.BULL_OPER_SA,sources.MODES,sources.ECJAN):
    #for fcsource in (sources.MODES,sources.ECJAN):
    #for fcsource in (sources.ECJAN,):
    for fcsource in (sources.BULL_OPER_SA,):

        checklastrun = False
        if checklastrun:
            print 'checking ', os.path.join(fcsource['dir'],fcsource['rexp'])
            a = sorted( glob.glob( os.path.join(fcsource['dir'],fcsource['rexp'])) )
            filelist[fcsource['shortname']],fcsource['lastrun'] = select.lastrun(a,date.year)

            dtg = fcsource['lastrun']
            date = dt.datetime.strptime(dtg,'%Y%m%d%H')

        ## EvdP keep date
        date = date 

        # now parse the fields extracted by Harmonie itself:
        # parse vfld files:
        cases = {'BULLSA':'36h14',
                 'MODES':'modes',
                 'ECJAN':''}
        cases_intar = {'BULLSA':'36h14',
                       'MODES':'modes',
                       'ECJAN':'ECJAN4'}

        tardir = fcsource['wdir'] #'/net/bhw284/nobackup/users/tijm/verg/HARM36'
        dbdir  = '/nobackup_1/users/plas/verif/BULL'
        #print date,tardir; sys.exit(1)

        if 0: # EvdP NB do not parse
            hb = 48
            # parses files (default 24 h) before date given:
            dbfile = parse_vfld.parse_dir(tardir,tempdir = tempdir,outdir = dbdir, #'/usr/people/plas',
                                          case    = cases[fcsource['shortname']],
                                          casetar = cases_intar[fcsource['shortname']],
                                          latest = False, date = date, hbefore = hb, hafter = 0,
                                          verbose = False) #NOT True for batch operation!
        else:
            dbfile = os.path.join(dbdir,'vfld_'+cases_intar[fcsource['shortname']]+'_'+dtg[0:6]+'.db')

        #sys.exit(1)


        # not exactly elegant: should be solved accordingly in websitegenerator.py
        outdir  = os.path.join(topdir,fcsource['shortname'],str(dtg[0:8]),str(dtg[8:10]),'time','0','field','nl')
        try:
            os.makedirs(outdir)
        except:
            print 'Outdir already there:',outdir

        for station in sid:

            ## outputfile
            outfile = os.path.join(outdir,'fc_ob_ts_'+fcsource['shortname']+'_'+dtg+'_'+str(station)+'.png')

            #extractdate = date - dt.timedelta(hours= 24); print 'extract date', extractdate
            ## extract modeldata from database
            extractdate = date; print 'extract date', extractdate,dbfile
            vdata = vis.extract_data(dbfile,extractdate,table='determ',s_id=station,hours = 36,hbefore = 36,hcycle = hc)
            for key in vdata.keys(): print 'leadtime,data ',key, len(vdata[key])
            #print vdata; sys.exit(0)

            #plotdate = date - dt.timedelta(hours= 24); print 'plotdate:',plotdate
            hb = 48 #36  # try making plots from hb hours before time it runs
            plotdate = date; print 'plotdate:',plotdate
            pf    = vis.create_plot(plotdate,vdata,hh,hb,hc,plotfile = outfile,s_id = station,
                                    synoplist = '/nobackup_1/users/plas/python/fabriek/fabriek/tools/allsynop.list')

        # switch for rest of stat figures
        continue

        ## include some HARP stats
        if fcsource['shortname'] == 'BULLSA':
            #for harps in ('BULL','BULL_NL'):
            #    if
            outdir  = os.path.join(topdir,fcsource['shortname'],str(dtg[0:8]),str(dtg[8:10]),'stats','0','field','nl')
            mkdir_p(outdir)
            figpath   = '/nobackup/users/plas/harp/BULL_NL/vis_harp_box*BULL_NL*'
            figpathnl = '/nobackup/users/plas/harp/BULL_3H/vis_harp_box*BULL_3H*'
            for sf in glob.glob(figpathnl):
                destsf = os.path.join(outdir,os.path.split(sf)[1])
                try:
                    os.symlink(sf, destsf)
                except OSError:
                    print 'link already exists, ',destsf

            outdir  = os.path.join(topdir,fcsource['shortname'],str(dtg[0:8]),str(dtg[8:10]),'stats','0','field','eur')
            mkdir_p(outdir)
            figpatheur = '/nobackup/users/plas/harp/BULL/vis_harp_box*BULL*'
            for sf in glob.glob(figpatheur):
                destsf = os.path.join(outdir,os.path.split(sf)[1])
                try:
                    os.symlink(sf, destsf)
                except OSError:
                    print 'link already exists, ',destsf
