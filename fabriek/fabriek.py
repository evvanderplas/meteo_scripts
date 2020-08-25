#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle
import time
import datetime as dt
import numpy
#import matplotlib, matplotlib.pyplot
#from mpl_toolkits.basemap import Basemap # map functionality
#from matplotlib import cm as CM

import pygrib

from tools.memory_inspector import report_mem_usage
backend = 'TkAgg'
backend = 'QT4Agg'
backend = 'Agg'
import matplotlib
matplotlib.use(backend)
import matplotlib.pyplot

#import h5tools,gribtools,finddate

def suf2(h):
    if h < 10:      return '0'+str(h)
    elif 10<=h<100: return str(h)
    

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

    import tools.read_stations as read_stations

    t0 = time.clock()
    today = dt.datetime.today()
    td = dt.datetime.strftime(today,'%Y%m%d%H')

    # all available maps at this point: 
    # nl:  netherlands
    # eur: large part of europe (BULL 800x800 domain)
    # rot: area around rotterdam
    # ijs: area around IJsselmeer

    #map_types   = ["nl", "eur", "rot", "ijs"] 
    #map_types   = ["nl", "eur", "rot"]
    map_types   = ["eur","nl"]

    # import which files to use for the fabriek from the file sources.py
    # e.g. ECJAN, ANJAN, BULL, Hirlam, etc
    import sources
    topdir = '/nobackup_1/users/plas/fabriek/'

    #debug:
    for source in sources.sources_list:
        print 'Using ',source['name'],map_types
    #print sources.obs_list #[0]['name']


    # make plots of the current hdf5 radar files:
    #plot_radar = True 
    #plot_radar = False

    # generate a list of (grib)files for the specified sources
    filelist = {}
    take_last_run = True

    for source in sources.sources_list:
        print source['shortname']
        
        # make a list of available files (must have at least +012 lt)
        sourcepath = os.path.join(source['dir'],source['rexp']); print sourcepath
        a = sorted( glob.glob( sourcepath) )
        filelist[source['shortname']] = a
        #for f in filelist[source['shortname']]: print f 
        #sys.exit(0)

        latestrun=0 # create variable that tracks latest run for obs directory
        if take_last_run:
            import tools.select as select
            
            filelist[source['shortname']],source['lastrun'] = select.lastrun(a,today.year)

            if int(source['lastrun']) > latestrun: latestrun = source['lastrun']

        print 'Latest run: ',latestrun; #sys.exit(0)

    if 0: # debug: which files?
        for source in sources.sources_list:
            print 'Checking: ',source['name'],source['lastrun']
            for grib_file in filelist[source['shortname']]: print grib_file
        sys.exit(1)


    # setting up parameter list to plot
    import param_lists,plottypes,settings

    TS = False # extract timeseries data
    if TS:
        import tools.timeseries as ts

        stationfile = 'tools/stations.txt'
        ts_dict = dict()
        tsfile  = './seriesdata.pkl' 

    # tool to check how much memory is consumed
    print 'starting memory management'
    report_mem_usage() # init

    for source in sources.obs_list:
        #plottypes.hdf_plot(source, map_types, latestrun, outdir = topdir)
        pass

    #sys.exit(1)

    for source in sources.sources_list:

        # could be file-source dependent (?)
        nllist = [sources.BULL_OPER_SA,sources.MODES,sources.RUC3dv,sources.fourdv]
        if source in nllist: 
            map_types = ['nl'] #,'rot']
        else:
            map_types = ['eur','nl']
            #pass
            

        if TS: # timeseries at stations:
            stations    = read_stations.stationdict(stationfile,filelist[source['shortname']][0], model = source)

            ts_dict[source['shortname']] = ts.stationlist()
            for station in stations:
                stat = stations[station]
                ts_dict[source['shortname']].add_station(stat['name'],stat['latlon'],stat['stationid'])

        grbindx_prev = None
        for grib_file in filelist[source['shortname']]: 
            
            print '************ open file ',grib_file
            print 'processing ',source['name'],time.clock() - t0, int(time.clock() - t0)*'='

            grbindx=pygrib.index(grib_file,
                                 'indicatorOfParameter',
                                 'indicatorOfTypeOfLevel',
                                 'level',
                                 'timeRangeIndicator')
            
            if TS:
                plottypes.point_extract(grbindx, source, stations, ts_dict[source['shortname']], outdir = topdir)
            #sys.exit(1)

            print 'Current gribs: ',grbindx, 'previous gribs: ',grbindx_prev

            plottypes.prec_plot(grbindx, grbindx_prev, source, map_types, outdir = topdir)
            print 'prec_plot: ',
            report_mem_usage() # report

            #grbindx_prev = grbindx; continue # debug..

            plottypes.windplot(grbindx, source, map_types, outdir = topdir)
            print 'windplot: ',
            report_mem_usage() # report

            plottypes.temp_plot(grbindx, source, map_types, outdir = topdir)
            print 'temp_plot: ',
            report_mem_usage() # report

            plottypes.cloud_plot(grbindx, source, map_types, outdir = topdir)
            print 'cloud_plot: ',
            report_mem_usage() # report

            plottypes.vis_plot(grbindx, source, map_types, outdir = topdir)
            print 'cloud_plot: ',
            report_mem_usage() # report

            # EvdP NB!
            #sys.exit(1) # one file done

            if grbindx_prev is not None:
                grbindx_prev.close()    # close previous(?)
                del(grbindx_prev)       # and deallocate
                
            grbindx_prev = grbindx  # def new previous


        #for k,v in globals().items(): print k,'=',v

        #sys.exit(1) # all files done

    if TS:
        f = open(tsfile,'wb')
        pickle.dump(ts_dict,f)
        f.close()

    sys.exit(1) # all sources done

