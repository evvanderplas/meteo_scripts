#! /usr/bin/env python

import os,sys
import datetime as dt
import nc_class_station as ncprof
import suites_info
import daterange

if len(sys.argv) < 5:
    print 'Usage: model, startingdate, indir, outdir'

    
    #sys.exit(1)
else:
    print sys.argv[1:6]
    model,stdate,indir,outdir = sys.argv[1:5]

startdate = dt.datetime(2010,4,2,0)
enddate   = dt.datetime(2013,1,1,0)

for d in daterange.daterange(startdate, enddate, delta = ('hours',12)):

    model,stdate,indir,outdir = 'REFO', ymdh, '/nobackup_1/users/plas/ncdat/prof/REFO/', '/nobackup_1/users/plas/ncdat/REFO'

    info   = suites_info.__dict__[model]
    stdate = d #dt.datetime.strptime(stdate,"%Y%m%d%H")
    nobj   = ncprof.ncconv(stdate,info=info,indir = indir,outdir=outdir)
    result = nobj.init_info(force = force) # force overrides skipping when output is already there
    if result != 'Done':
        #nobj.loop_over_grib()
        nobj.initprofiles()
        print 'Writing files'
        nobj.writefile()
        print 'Post-processing'
        nobj.postproc_uv()
        nobj.postproc_pz()
        nobj.postproc_tq()
        nobj.postproc_radiation()
        nobj.postproc_percent()
        print 'Closing files'
        nobj.wrapup()

else:
    print 'Something wrong:'
    print 'Usage: model, startingdate, indir, outdir'
    print sys.argv[1:]
    print model, model in suites_info.__dict__.keys(), suites_info.__dict__.keys()
