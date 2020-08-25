#! /usr/bin/env python

import os,sys,glob
import numpy as np
import datetime as dt
import subprocess
#import tarfile
from os.path import join as pjoin

import netCDF4
import pygrib

#import nc_profile_class as ncprof
import nc_class as ncprof
import daterange
import untar
import suites_info_BHW477 as suites_info

suites  = suites_info.__dict__
mysuite = suites['BULL_REFO']
outdir  = '/net/bhw379/nobackup_1/users/plas/ncdat/prof'

if __name__ == '__main__':

    print sys.argv

    cli_cmd  = 'python /usr/people/plas/python/tools/use_ncprof_cli.py {model} {startdate} {indir} {outdir}'
    surf_cmd = '/usr/people/plas/python/tools/convert_harmonie_ext.py -n -s -l {levs} -t {tab} -prt -O -o {outfile} {inf}'

    ## create datetime of latest possible run:
    begindate = dt.datetime(2010,4,21,0)
    enddate   = dt.datetime(2010,9,1,0)

    # bhw477
    #begindate = dt.datetime(2011,3,5,0)
    #enddate   = dt.datetime(2012,1,1,0)
    #enddate   = dt.datetime(2011,3,4,12)
    #enddate   = dt.datetime(2011,7,31,12)

    #begindate = dt.datetime(2010,9,13,0)
    #enddate   = dt.datetime(2011,1,1,0)

    # bhw478
    begindate = dt.datetime(2011,10,3,0)
    enddate   = dt.datetime(2011,12,31,12)

    for sdate in daterange.daterange(begindate,enddate): #,delta=('hours',12)):

        for st in range(0,25,12):

            # today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0,microsecond=0)
            stime = sdate + dt.timedelta(hours=st)
            #stime = dt.datetime(2014,6,1,0) #NB EvdP debug!
            ymdh  = dt.datetime.strftime(stime,'%Y%m%d%H')
            print stime,ymdh

            # untar
            fctar = pjoin(mysuite['tardir'],mysuite['tgzname'].format(ymdh=ymdh))
            print 'Untarring ',fctar
            if not os.path.exists(fctar): 
                print 'Not there',fctar
                continue
                
            gdir,fformat,ltmax = mysuite['dir'],mysuite['fformat'],mysuite['ltmax']
            lastfile = pjoin(gdir,fformat.format(ymdh = ymdh,lt = str(ltmax).zfill(3)))
            if not os.path.exists(lastfile):
                untar.untar(fctar,mysuite['dir'],pattern=mysuite['fformat'].format(ymdh=ymdh,lt='*'),tempdir=pjoin(mysuite['dir'],'untar'))
            #sys.exit(0)

            # check if untarring went ok: otherwise skip to next
            if not os.path.exists(lastfile): continue

            # surface:
            inf     = pjoin(gdir,fformat.format(ymdh = ymdh,lt='*'))
            outfile = pjoin(outdir,mysuite['shortname'],'harm_REFO_{ymdh}.nc'.format(ymdh = ymdh))
            cmd = surf_cmd.format(tab=mysuite['tab'], outfile = outfile, inf = inf,levs=mysuite['nlevs'])
            print cmd

            # run command
            pid = subprocess.Popen(cmd,shell=True)
            pid.communicate()
            ##os.system(cmd)

            # profile:
            cmd = cli_cmd.format(model = 'BULL_REFO', 
                                 startdate=ymdh,
                                 indir = mysuite['dir'], 
                                 outdir = os.path.join(outdir,mysuite['shortname'])) 
            
            print cmd

            # run command
            pid = subprocess.Popen(cmd,shell=True)
            pid.communicate()
            

            for f in sorted(glob.glob(pjoin(mysuite['dir'],mysuite['fformat'].format(ymdh=ymdh,lt='*')))):
                print 'Remove',f, os.path.exists(f)
                os.remove(f)

