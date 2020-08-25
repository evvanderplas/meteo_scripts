#! /usr/bin/env python

import os,sys
import datetime as dt
import nc_class_station as ncprof
import suites_info
import argparse

if len(sys.argv) < 5:
    print 'Usage: model, startingdate, indir, outdir'

    model,stdate,indir,outdir = 'VETREFO', '2010071406', '/nobackup_1/users/plas/temp/2010/07/14/06/', './'
    #sys.exit(1)
else:
    print sys.argv[1:6]
    model,stdate,indir,outdir = sys.argv[1:5]

#if len(sys.argv) >= 6
force = False
force = True # overwrite, especially for debugging

if model in suites_info.__dict__.keys():

    info   = suites_info.__dict__[model]
    stdate = dt.datetime.strptime(stdate,"%Y%m%d%H")
    nobj   = ncprof.ncconv(stdate,info=info,indir = indir,outdir=outdir)
    result = nobj.init_info(force = force) # force overrides skipping when output is already there
    if result != 'Done':
        nobj.loop_over_grib()
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
