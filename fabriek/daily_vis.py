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
    dtg  = date.strftime('%Y%m%d%H')


    #for fcsource in (sources.BULL_OPER_SA,sources.MODES,sources.ECJAN): # EvdP ANJAN abandoned
    #for fcsource in (sources.BULL_OPER_SA,sources.MODES,sources.RUC3dv): # EvdP MODES temporarily (?) shut down
    for fcsource in (sources.BULL_OPER_SA,sources.fourdv,sources.BULL_38h11,):
    #for fcsource in (sources.BULL_38h11,):

        today = dt.datetime.today()
        today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0)
        for h in range(0,25,3):
            date = today - dt.timedelta(hours=h)
            y,m,d,h = date.strftime('%Y'),date.strftime('%m'),date.strftime('%d'),date.strftime('%H')
            
            print 'checking ', os.path.join(fcsource['dir'].format(y=y,m=m,d=d,h=h),fcsource['rexp']) 
            
        a = sorted( glob.glob( os.path.join(fcsource['dir'],fcsource['rexp'])) )
        filelist[fcsource['shortname']],fcsource['lastrun'] = select.lastrun(a,date.year)
        dtg  = fcsource['lastrun']
        print dtg
        try:
            date = dt.datetime.strptime(dtg,'%Y%m%d%H')

            #y,m,d,h = date.strftime('%Y'),date.strftime('%m'),date.strftime('%d'),date.strftime('%H')
            #print 'checking ', os.path.join(fcsource['dir'].format(y=y,m=m,d=d,h=h),fcsource['rexp'])
        except:
            print 'Date from lastrun not working'

        #print ; #sys.exit(1)
        ## EvdP try yesterday:
        #date = date 

        # now parse the fields extracted by Harmonie itself:
        # parse vfld files:
        cases = {'BULLSA':'36h14',
                 'MODES':'modes',
                 'ECJAN':'',
                 'RUC':'ruc',
                 '4DVAR':'4DVAR_300',
                 '38H11':'38h11',
                 }
        cases_intar = {'BULLSA':'36h14',
                       'MODES':'modes',
                       'ECJAN':'ECJAN4',
                       'RUC':'ruc',
                       '4DVAR':'4DVAR_300',
                       '38H11':'38h11',
                   }

        tardir = fcsource['wdir'] #'/net/bhw284/nobackup/users/tijm/verg/HARM36'
        #print date,tardir; sys.exit(1)

        if 1: # Update database (EvdP NB)
            hb = 48
            #hb = 144
            # parses files (default 24 h) before date given:
            dbfile = parse_vfld.parse_dir(tardir,tempdir = tempdir,outdir = '/nobackup_1/users/plas/verif/BULL', #'/usr/people/plas',
                                          case    = cases[fcsource['shortname']],
                                          casetar = cases_intar[fcsource['shortname']],
                                          latest = False, date = date, hbefore = hb, hafter = 0,
                                          verbose = False) #NOT True for batch operation!
        else:
            # do not update database
            dbdir = '/nobackup_1/users/plas/verif/BULL'
            dbfile = os.path.join(dbdir,'vfld_'+cases_intar[fcsource['shortname']]+'_'+dtg[0:6]+'.db')
            #print 'Ja?', dbfile,os.path.exists(dbfile); sys.exit(1)

        #sys.exit(1)


        # not exactly elegant: should be solved accordingly in websitegenerator.py
        outdir  = os.path.join(topdir,fcsource['shortname'],str(dtg[0:8]),str(dtg[8:10]),'time','0','field','nl')

        # first attempt:
        outdir  = os.path.join(topdir,fcsource['shortname'],'time')
        archdir = os.path.join(topdir,fcsource['shortname'],'time','old')

        try:
            os.makedirs(outdir)
        except:
            print 'Outdir already there:',outdir

        try:
            os.makedirs(archdir)
        except:
            print 'archive already in place',archdir
            
        if 1: #try:
            for f in glob.glob(os.path.join(outdir,'*png')):
                print f
                os.rename(f,os.path.join(archdir,os.path.split(f)[1]))
        else: #except:
            print 'no previous images? check ',outdir

        #sys.exit(0)

        for station in sid:

            ## outputfile
            outfile = os.path.join(outdir,'fc_ob_ts_'+fcsource['shortname']+'_'+dtg+'_'+str(station)+'.png')

            #extractdate = date - dt.timedelta(hours= 24); print 'extract date', extractdate
            ## extract modeldata from database
            extractdate = date; print 'extract date', extractdate,dbfile
            vdata = vis.extract_data(dbfile,extractdate,table='determ',s_id=station,hours = 36,hbefore = 36,hcycle = hc)
            for key in vdata.keys(): print key, len(vdata[key])

            #plotdate = date - dt.timedelta(hours= 24); print 'plotdate:',plotdate
            hb = 48 #36  # try making plots from hb hours before time it runs
            plotdate = date; print 'plotdate:',plotdate
            pf    = vis.create_plot(plotdate,vdata,hh,hb,hc,plotfile = outfile,s_id = station,
                                    synoplist = '/nobackup_1/users/plas/python/fabriek/fabriek/tools/allsynop.list')


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
