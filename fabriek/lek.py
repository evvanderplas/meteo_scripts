#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle
import datetime as dt
import numpy


import pygrib

backend = 'TkAgg'
backend = 'QT4Agg'
import matplotlib
matplotlib.use(backend)
import matplotlib.pyplot

import resource


###############################################
#
# Main program
#
###############################################

if __name__ == '__main__':

    map_types   = ["nl", "eur", "rot"]
    map_types   = ["eur"]

    # import which files to use for the fabriek
    # e.g. ECJAN, ANJAN, BULL, Hirlam, etc
    import sources
    for source in sources.sources_list:
        print 'Using ',source['name'],
 
    filelist = []
    for source in sources.sources_list:
        #print source
        a = sorted( glob.glob( os.path.join(source['dir'],source['rexp'])) )
        filelist.append( a )


    # setting up parameter list to plot
    import param_lists,settings,plottypes

    #param_list,    searchlist   = param_lists.genSearchList(settings.parameter_list)
    postparam_list,ppsearchlist = param_lists.genSearchList(settings.postpr_list)

    # debug memory leakage
    import gc
    gc.enable()
    #gc.set_debug(gc.DEBUG_LEAK)

    #print searchlist;sys.exit(1)
    for s,source in enumerate(sources.sources_list):

        # could be file-source dependent (?)
        outdir = './'
        if source['shortname'] == 'EXP': outdir = source['outdir']
        #print outdir; sys.exit(1)

        grbindx_prev = None
        for grib_file in filelist[s][:5]:
            
            print '************ open file ',grib_file

            grbindx=pygrib.index(grib_file,
                                 'indicatorOfParameter',
                                 'indicatorOfTypeOfLevel',
                                 'level')
            
            #plottypes.windplot(grbindx, source, map_types, outdir = outdir)

            plottypes.temp_leak_plot(grbindx, source, map_types, outdir = outdir)
            print 'SELF ',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print 'CHILDREN', resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

            #plottypes.cloud_plot(grbindx, source, map_types, outdir = outdir)

            #plottypes.prec_plot(grbindx, grbindx_prev, source, map_types, outdir = outdir)
            #sys.exit(1) # one file done

            if grbindx_prev is not None: grbindx_prev.close()    # close previous(?)
            grbindx_prev = grbindx  # def new previous
            #grbindx.close()

             
            #for k,v in globals().items(): print 'globals:', k,'=',v

        #sys.exit(1) # all files done

    #sys.exit(1) # all sources done
    gc.collect()

    import objgraph
    objgraph.get_leaking_objects()
