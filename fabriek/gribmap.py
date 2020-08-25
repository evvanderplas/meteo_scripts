#!/usr/bin/env python
import argparse, os, sys, glob
import pygrib

#import fabriek.sources
#import fabriek.plottypes as fabplots
import os,sys,re
import sources, settings
import get_plotinfo as getinfo
import plottypes 

# characters not allowed in strings for file names:
re_ill =  re.compile(r"^[^<>/{}[\]~`\s\*\\\?\|]*$");

if __name__ == '__main__':

    print sys.argv
    parser = argparse.ArgumentParser(prog = 'GRIBMAP',description='Parse plotting arguments.')
    parser.add_argument('-v', '--verbose',
                        default=False, action='store_true',
                        help="Give more verbose output." )
    parser.add_argument('-o', '--outdir',  
                        type=str,   default='./',
                        action='store', 
                        dest='outdir', 
                        #required=True, 
                        help="The output directory" )
    parser.add_argument('-p', '--parameter',    
                        type=int,   default=None, #61,    
                        action='store', 
                        dest='param',
                        help='The plottable grib parameter (e.g. 11 for temperature,61 for precip, depends on table)')
    parser.add_argument('-lt', '--leveltype',    
                        type=int,   default=None, #456,    
                        action='store', 
                        dest='leveltype',
                        help="The leveltype: 105 for surface level, 109 for hybrid model level, 103 for sigma level, etc (kpds6)" )    
    parser.add_argument('-l', '--level',    
                        type=int,   default=None, #456,    
                        action='store', 
                        dest='level',
                        help="Vertical model level" )    
    parser.add_argument('-tr', '--timerange',    
                        type=int,   default=None, #0 or 4,    
                        action='store', 
                        dest='tr_ind',
                        help="Time range indicator (accumulated:4 or intensity:0  etc)" )    
    parser.add_argument('-s', '--setting',    
                        type=str,   default= None, #'APCP',    
                        action='store', 
                        dest='sett',
                        help="Type of plottable parameter (APCP,T2M,PRES, etc)" )    
    parser.add_argument('-t', '--tab',    
                        type=str,   default= None, # Harmonie 36 or 37?
                        action='store', 
                        dest='gtab',
                        help="Grib table, eg har36, har37, hir, ..." )    
    parser.add_argument('-d', '--domain',    
                        type=str,   default='nl',    
                        action='store', 
                        dest='domain',
                        help="Domain on which to plot: for now nl (Netherlands),eur (Europe),rot (Rotterdam)")
    parser.add_argument('-c', '--case',    
                        type=str,   default='HARM',    
                        action='store', 
                        dest='case',
                        help="name of case or model to put in name of figure etc to distinguish between different series of plots"
                        )
    parser.add_argument('infiles', metavar='inputfiles', 
                        type=str,      
                        action='append',
                        nargs='*',
                        help='(GRIB) input files to plot')

    args = parser.parse_args()

    # check input file(s)
    print args.infiles #, not args.infiles
    if args.infiles == [[]]:
        parser.print_help()
        assert False, "no input files"


    # see if it is necessary to force the parameter to what you want it to be.
    force = False
    if args.sett and args.param: 
        force = True

    # allow for sloppy settings:
    if args.sett:
        if args.sett in ('PCP','APCP','NCPCP','precip','RR'):
            args.sett = settings.PCP_acc
        elif args.sett.lower() in ('pcpi',):
            args.sett = settings.PCP_i
        elif args.sett.lower() in ('vv','w',):
            args.sett = settings.VV
        elif args.sett.lower() in ('t2m','t','temp2m','temp','temperature'):
            args.sett = settings.T2M
        elif args.sett in ('U10','V10','uv','UV','wind'):
            args.sett = settings.WIND
        elif args.sett.lower() in ('gust',):
            args.sett = settings.GUST
        elif args.sett in ('P','MSLP','pres','PRES','press','PRESS','pressure'):
            args.sett = settings.PRESS
        elif args.sett in ('CLC','clc','cloud','CLOUD','cloudcover'):
            args.sett = settings.CLOUD
        elif args.sett in ('VIS','vis','mist','soep','zicht'):
            args.sett = settings.VIS
    elif args.sett == None:
            args.sett = settings.default

    print 'setting:',args.sett

    if args.param:
        # check if there is a setting with that grib parameter:
        for s in settings.all_list:
            if args.level is not None:
                args.sett['level'] = args.level
                if args.param == s['grib_indicator'] and args.sett['level'] == s['level']:
                    print 'processing as if ',s
                    yn = raw_input('processing as if '+s['name']+' (Y/n)')
                    if yn in ('y','Y','yes','ok','\n',''):
                        args.sett = s
                        break
                    else:
                        print 'Ni';sys.exit(1)
                
            if args.param == s['grib_indicator']:
                print 'processing as if ',s
                yn = raw_input('processing as if '+s['name']+' (Y/n)')
                if yn in ('y','Y','yes','ok','\n',''):
                    args.sett = s
                    break
                elif yn == 'n':
                    pass
                else:
                    print 'Ni';sys.exit(1)

        args.sett['grib_indicator'] = args.param

    if args.param is not None:
        args.sett['grib_indicator'] = args.param
        print 'Setting ',args.sett['name'], ' with param ',args.sett['grib_indicator']

    if args.level is not None:
        args.sett['level'] = args.level

    if args.leveltype is not None and type(args.leveltype) == int: 
        args.sett['level_indicator'] = args.leveltype
        print 'YES: ',args.leveltype, args.sett['level_indicator']

    if args.tr_ind is not None:
        args.sett['tr_ind'] = args.tr_ind

    if args.sett not in settings.all_list:
        #for s in  settings.all_list: print s
        print 'Setting not known: ',args.sett
        sys.exit(1)


    if args.domain not in ('nl','eur','ijs','rot'):
        print 'Domain not known: ',args.domain,', choose nl,eur,rot or ijs'
        sys.exit(1)

    # check if name is well-behaved:
    if re_ill.match(args.case):
        print("RE2: All chars are valid.",args.case)
        pass
    else:
        print("RE2: Not all chars are valid.",args.case)
        sys.exit(1)


    if args.verbose:
        print "Input:              ", args.infiles
        print "Output:             ", args.outdir #file
        print "parameter           ", args.param
        print "level               ", args.level
        print "time range          ", args.tr_ind
        print "table               ", args.gtab
        print "setting             ", args.sett
        print "domain              ", args.domain
        print "case                ", args.case

    flattend = []
    if args.infiles:
        for sublist in args.infiles:
            flattend += sublist 
        args.infiles = flattend

    #print args.sett
    #sys.exit(1)

    if os.path.splitext(args.infiles[0])[1] == '.h5':
        print os.path.split(args.infiles[0])[1][:8]
        if os.path.split(args.infiles[0])[1][:8] == 'RAD_NL23':
            source = sources.RADHDF23
        elif os.path.split(args.infiles[0])[1][:12] == 'RAD_NL25_RAC':
            source = sources.RADHDF25_3H
        elif os.path.split(args.infiles[0])[1][:8] == 'RAD_NL25':
            source = sources.RADHDF25
        elif 'CMa' in os.path.split(args.infiles[0])[1]:
            source = sources.MSG
        else:
            print 'Assuming it is dutch radar...'
            source = sources.RADHDF25
    else:
        if args.gtab == 'hir':
            source  = sources.HIRLAM_D11            
        elif args.gtab == 'har37':
            source  = sources.BULL_OPER
            source['version'] = 3712
        else: # args.gtab == 'har37':
            source  = sources.BULL_OPER
            source['version'] = 3614

    #domain  = 'nl' #'eur'
    #setting = fabriek.settings.PCP_acc
    #setting = settings.PCP_acc ; #print setting # args.sett

    if 1:

        import numpy 

        print args.infiles #os.path.splitext(args.infiles[0])
        #if os.path.splitext(args.infiles[0])[1] == '.h5': # args.infiles[0][-2:] == 'h5':
        if source in (sources.RADHDF23,sources.RADHDF25,sources.MSG): 
            for h5file in args.infiles:
                h5info = getinfo.get_hdf_plotinfo(h5file,args.sett,source)

                # here I concatenate the "case" string to the existing
                plottypes.contourplot(h5info,domain = args.domain,outdir = args.outdir,lsmask_colour = 'red')

        else:
            
            grbindx_prev = None
            for grib_file in args.infiles:
                
                grbindx=pygrib.index(grib_file,
                                     'indicatorOfParameter',
                                     'indicatorOfTypeOfLevel',
                                     'level',
                                     'timeRangeIndicator')

                print 'before:',args.sett
                var   = getinfo.get_grib_plotinfo(grbindx,args.sett,source,force = force)
                print 'After getting info: ',var
                #sys.exit(1)

                if args.sett == settings.VIS:

                    plottypes.vis_plot(grbindx,source,[args.domain],outdir = args.outdir,ftree = False)
                    continue

                    #getinfo.prec_plot(grbindx,grbindx_prev,source,[args.domain],outdir = args.outdir)
                    tpcp = getinfo.get_grib_plotinfo(grbindx,settings.MLRA,source)
                    tclw = getinfo.get_grib_plotinfo(grbindx,settings.CLWAT,source)
                    tgrp = getinfo.get_grib_plotinfo(grbindx,settings.MLGR,source)
                    tsnw = getinfo.get_grib_plotinfo(grbindx,settings.MLSN,source)

                    dl = 0.
                    #dl = 1.e-10
                    extinction = (
                        144.7 * numpy.power(numpy.select([tclw['values']<0.,tclw['values']>=0],[dl,1.2e3*tclw['values']]  ),0.88) + 
                        1.1   * numpy.power(numpy.select([tpcp['values']<0.,tpcp['values']>=0],[dl,1.2e3*tpcp['values']]),0.75) +
                        10.4  * numpy.power(numpy.select([tsnw['values']<0.,tsnw['values']>=0],[dl,1.2e3*tsnw['values']]),0.78) + 
                        2.4   * numpy.power(numpy.select([tgrp['values']<0.,tgrp['values']>=0],[dl,1.2e3*tgrp['values']]),0.78)
                        )

                    var['values'] = 1000.*3.912 / (extinction + dl)

                elif args.sett == settings.WIND:
                    u = getinfo.get_grib_plotinfo(grbindx,settings.U10,source)
                    v = getinfo.get_grib_plotinfo(grbindx,settings.V10,source)
                    var['values'] = numpy.sqrt(u['values']**2 + v['values']**2)

                elif args.sett == settings.GUST:
                    u = getinfo.get_grib_plotinfo(grbindx,settings.UGST,source)
                    v = getinfo.get_grib_plotinfo(grbindx,settings.VGST,source)
                    var['values'] = numpy.sqrt(u['values']**2 + v['values']**2)
                    

                #print var['values']
                #for item in accprec.keys(): print item, accprec[item]

                if var is not None:

                    var['source_short'] = args.case

                    ## functionality should be implemented elsewhere, I suppose...
                    if var['param'] == 11 and var['values'].max() > 200: # convert Kelvin to Celsius
                        var['values'] -= 273.15
                    #elif var['param'] in (33,34): # convert m/s to knots
                    #    var['values'] *= 1.9438444924406046
                    elif var['param'] in (61,62,63,181,184,201) and var['level'] in (0,457): # convert to mm/hour
                        if (var['param'] in (181,184,201) and var['level'] in (0,) and var['tr_ind'] == 0): # new Harmonietab
                            print 'precip intensity'
                            var['values'] *= 3.6e3
                        else:
                            try:
                                var['values'] -= var_prev
                                print 'deducted previous timestep'
                            except NameError: 
                                print 'no values to deduct'
                                pass
                            var_prev = var['values']
                    elif (var['param'] == 61 and var['level'] == 456):
                        print 'precip intensity'
                        var['values'] *= 3.6e3

                    elif var['param'] in (142,143): # ECMWF
                        var['values'] *= 1.e3


                    plottypes.contourplot(var,domain = args.domain,outdir = args.outdir)

                else:
                    print 'No valid plotting information found'

                #grbindx_prev = grbindx
                # end of for loop over files
