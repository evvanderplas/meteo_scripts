#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle, shutil
import time
import datetime as dt
import numpy
#import matplotlib, matplotlib.pyplot
#from mpl_toolkits.basemap import Basemap # map functionality
#from matplotlib import cm as CM
sys.path.append('/usr/people/plas/python/tools')

import pygrib

from tools.memory_inspector import report_mem_usage
backend = 'TkAgg'
backend = 'QT4Agg'
backend = 'Agg'
import matplotlib
matplotlib.use(backend)
import matplotlib.pyplot

#import h5tools,gribtools,finddate

def hour_range(begindate,leadtime = 48,dl=1,enddate = None):

    if type(begindate) != dt.datetime: 
        print 'hour_range should be called with datetime object, got',begindate
        sys.exit(1)

    begindate = begindate.replace(minute=0,second=0,microsecond=0)

    for h in range(0,leadtime,dl):
        yield begindate + dt.timedelta(hours=h)

def clean_obsdir(outdir):

    for d in  os.listdir(os.path.join(outdir,'obs')):
        #print d
        if d[0:3] == '201': # it is an imagedirectory with softlinks, DELETE it
            print 'Deleting old softlinks: rm -rf {d}'.format(d=d)
            shutil.rmtree(os.path.join(outdir,'obs',d))

    return 0

###############################################
#
# Main program
#
# fabriek.py is a script that creates (png) plots of any available grib files
# * collects necessary files (specified in sources.py)
# * looks in settings.py which variables should be plotted, using what settings
# * * in plottypes.py are the scripts that extract necessary infor from grib file
#     and combine them into the plottable quantity
# * plots precipitation radar files
# * 
#
###############################################

#for file in gribec_files: print os.path.split(file)[1]

#settingsfile = './fabriek_settings.txt' 

if __name__ == '__main__':

    import bmap
    import tools.read_stations as read_stations
    sys.path.append('/usr/people/plas/python/tools')
    import plotinfo

    # import which files to use for the fabriek from the file sources.py
    # e.g. ECJAN, ANJAN, BULL, Hirlam, etc
    import sources

    single_exp,only_obs = False,False
    if len(sys.argv) > 1:
        if sys.argv[1] == 'obs': only_obs = True
        else: only_obs = False

        #print sources.__dict__.keys()
        if sys.argv[1] in sources.__dict__.keys():
            sources.sources_list = (sources.__dict__[sys.argv[1]],)
            single_exp = True
            print 'Processing only experiment ',sys.argv[1]
        if only_obs: print 'Processing only observations',only_obs,sources.obs_list
        #sys.exit(1)

    t0 = time.clock()

    ## create datetime of latest possible run:
    today = dt.datetime.today()
    today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0,microsecond=0)
    td = dt.datetime.strftime(today,'%Y%m%d%H')

    # all available maps at this point: 
    # nl:  netherlands
    # eur: large part of europe (BULL 800x800 domain)
    # rot: area around rotterdam
    # ijs: area around IJsselmeer

    #map_types   = ["nl", "eur", "rot", "ijs"] 
    #map_types   = ["nl", "eur", "rot"]
    map_types   = ["eur","nl"]

    topdir = '/nobackup_1/users/plas/fabriek/'

    #debug:
    for source in sources.sources_list:
        print 'Using ',source['name'],map_types
    print 'Obs: ',sources.obs_list #[0]['name']


    # make plots of the current hdf5 radar files:
    plot_radar = True 
    plot_radar = False

    # generate a list of (grib)files for the specified sources
    filelist = {}
    take_last_run = True

    # setting up parameter list to plot
    import param_lists,settings
    import plottypes_dev as plottypes

    # tool to check how much memory is consumed
    print 'starting memory management'
    report_mem_usage() # init


    clean_obsdir(topdir)

    if not single_exp:
        for source in sources.obs_list:

            if source == sources.RADAR:
                hdfmap_types = ['nl',]
            else:
                hdfmap_types = ['nl','eur']
            plottypes.hdf_plot(source, hdfmap_types, today, outdir = topdir)

        if only_obs: sys.exit(1)

    for source in sources.sources_list:

        # could be file-source dependent (?)
        if source in (sources.BULL_OPER_SA,sources.fourdv,sources.BULL_38h11):
            map_types = ['nl'] 
        elif source in (sources.BULL_OPER,):
            map_types = ['eur'] 
        else:
            map_types = ['eur','nl']
            #pass

        bdate = None
        print today
        for delay in range(0,24,3):
            st = today - dt.timedelta(hours=delay)
            gribdir  = source['dir'].format(y=st.strftime('%Y'),
                                           m=st.strftime('%m'),
                                           d=st.strftime('%d'),
                                           h=st.strftime('%H'))
            filename = source['fformat'].format(dtgh=st.strftime('%Y%m%d%H'),lt = '012')
            ff = os.path.join(gribdir,filename)
            
            print ff, os.path.exists(ff)
            if os.path.exists(ff):
                bdate = st
                break
        if bdate == None: 
            print 'No files of last 24 hours for ',source['name']
            break
            


        grbindx_prev = None
        for h in range(0,49,1): #hour_range(bdate):

            # construct gribfile from "dir" and "fformat" in sources:
            # source['dir'] is the directory where the files live, and 
            # source['fformat'] is the name of the file where yyyymmddhh 
            # and the leadtime are to be filled in
            #
            filename  = source['fformat'].format(dtgh=st.strftime('%Y%m%d%H'),lt = str(h).zfill(3))
            grib_file = os.path.join(gribdir,filename)
            
            print '************ open file ',grib_file, os.path.exists(grib_file)

            # if the file with leadtime h does not exist, jump to next source
            # no use to look further for files, I suppose
            if not os.path.exists(grib_file): break 

            print 'processing ',source['name'],time.clock() - t0, int(time.clock() - t0)*'='

            # open grib file, make index:
            grbindx=pygrib.index(grib_file,
                                 'indicatorOfParameter',
                                 'indicatorOfTypeOfLevel',
                                 'level',
                                 'timeRangeIndicator')
            
            # do the plotting: in plottypes are the "recipes" to make the necessary plots

            plottypes.prec_plot(grbindx, grbindx_prev, source, map_types, outdir = topdir)
            print 'prec_plot: ',
            report_mem_usage() # report

            plottypes.windplot(grbindx, source, map_types, outdir = topdir)
            print 'windplot: ',
            report_mem_usage() # report

            plottypes.temp_plot(grbindx, source, map_types, outdir = topdir)
            print 'temp_plot: ',
            report_mem_usage() # report

            brtemp = plottypes.cloud_plot(grbindx, source, map_types, outdir = topdir)
            print 'cloud_plot: ',
            report_mem_usage() # report

            plottypes.vis_plot(grbindx, source, map_types, outdir = topdir)
            print 'cloud_plot: ',
            report_mem_usage() # report


            #sys.exit(1) # one file done

            if grbindx_prev is not None:
                grbindx_prev.close()    # close previous(?)
                del(grbindx_prev)       # and deallocate
                
            grbindx_prev = grbindx  # define new previous gribindex, for accumulations and such


        #sys.exit(1) # all files done

    sys.exit(1) # all sources done

