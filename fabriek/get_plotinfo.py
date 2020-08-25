#! /usr/bin/env python

import os,sys,numpy
import datetime as dt

import param_lists


def init_list(model = 'Harmonie', version = '3712'):

    '''
    for different versions of Harmonie, different parameter codes are used:
    so either read a new table, or get some parameters for the old version
    new: use 181,184 etc for rain, snow etc, use timeRangeIndicator to discern between acc/instantaneous
    old: use 61,62 etc for rain,snow etc, but different for upper air, acc = level 457, inst = 456...
    Hirlam: ?
    '''

    scriptlocation =  os.path.split(os.path.realpath(sys.argv[0]))[0] #split(sys.argv[0])[0]
    

    print model,version

    gribtable = {}
    if model == 'Harmonie' and int(version) >= 3711: 
        tabfile = os.path.join(scriptlocation,'tools/2.233.253.table')
        #print os.getcwd(), scriptlocation, os.path.exists(tabfile); sys.exit(1)
        #tabfile = 'tools/2.233.253.table'

    elif model == 'Hirlam':
        tabfile = os.path.join(scriptlocation,'tools','2.99.hir.table')

    else:
        # use old gribtab
        tabfile = os.path.join(scriptlocation,'tools','2.99.128.table')

    #print version
    #print os.getcwd(), scriptlocation, os.path.exists(tabfile); sys.exit(1)

        #return gribtable

    f= open(tabfile,'r')
    for line in f:
        s = line.split()
        gribtable[s[1]] = int(s[0])
    f.close()

    gribtable['cloudlevel'] = 0

    gribtable['gustlevel'] = 10
    gribtable['gustlt']    = 105
    gribtable['gusttr'] = 0

    # not elegant: set (level,timeRangeIndicator)
    if model == 'Harmonie' and int(version) < 3711 :
        gribtable['acclevel'] = 457
        gribtable['intlevel'] = 456

        # timeRangeIndicator not used in Harmonie < 37h1
        gribtable['inttr']  = 0
        gribtable['acctr']  = 0

    elif model == 'Hirlam': # no gust output from Hirlam
        gribtable['acclevel'] = 0
        gribtable['intlevel'] = 0

        gribtable['acctr']  = 4
        gribtable['inttr']  = 0

    else:
        gribtable['acclevel'] = 0
        gribtable['intlevel'] = 0

        
        gribtable['acctr']  = 4
        gribtable['inttr']  = 0
        gribtable['gusttr'] = 2

    return gribtable

def complete_setting(setting,source,force=False):

    '''
    based on the source modify the entries as put in settings.py
    force = True makes it keep the settings as entered

    EvdP Awful solution, redesign...
    '''

    if source.has_key('version'):
        version =  source['version']
    else:
        version = 0
    gribtable = init_list(model = source['model'],version = version)

    if setting['time'] == 'ml': # for now just for lowest model levels:
        # due to wrong use of grib table, historically:
        if int(version) < 3711 :
            gribtable['r']  = 62
            gribtable['gr'] = 201
            gribtable['s']  = 79

        lt = 109
        gribtable['acclevel'] = 60
        gribtable['acctr']    = 0
    else:
        lt = 105

    replacedict = {
        'pcpi':   {'grib_indicator':gribtable['r'],   'level_indicator':lt,'level':gribtable['intlevel'],'tr_ind':gribtable['inttr']},
        'apcp':   {'grib_indicator':gribtable['r'],   'level_indicator':lt,'level':gribtable['acclevel'],'tr_ind':gribtable['acctr']},
        'snow':   {'grib_indicator':gribtable['s'],   'level_indicator':lt,'level':gribtable['acclevel'],'tr_ind':gribtable['acctr']},
        'graup':  {'grib_indicator':gribtable['gr'],  'level_indicator':lt,'level':gribtable['acclevel'],'tr_ind':gribtable['acctr']},
        'ugst'  : {'grib_indicator':gribtable['162'], 'level_indicator':lt,'level':gribtable['gustlevel'],'tr_ind':gribtable['gusttr']},
        'vgst'  : {'grib_indicator':gribtable['163'], 'level_indicator':lt,'level':gribtable['gustlevel'],'tr_ind':gribtable['gusttr']},
        'tcl': {'grib_indicator':gribtable['114'], 'level_indicator':lt,'level':gribtable['cloudlevel'],'tr_ind':gribtable['inttr']}
        }

    if force:
        if setting.has_key('grib_indicator'):
            if setting.has_key('level_indicator'):
                if setting.has_key('level'):
                    return setting
                else:
                    print 'solve later: forced: ',setting; sys.exit(1)

    elif setting['shortname'] in replacedict.keys() and int(version) > 3711:
        print 'replacing', setting['shortname'],replacedict[setting['shortname']]
        for k,v in replacedict[setting['shortname']].iteritems():
            setting[k] = v
        if setting['time'] == 'int':
            setting['level']  = gribtable['intlevel']
            setting['tr_ind'] = gribtable['inttr']
        elif setting['time'] == 'ml':
            setting['level']      = 60 
            setting['level_indicator']  = 109
            setting['tr_ind']     = 0

    return setting

def get_grib_plotinfo(grib_index,settings,source,par = 'acc',force = False):


    print 'In get_grib_plotinfo: first ',settings['grib_indicator'],settings['level_indicator'],settings['level']

    #from settings import init_list
    settings['time'] = par
    settings = complete_setting(settings,source,force = force)


    param     = settings['grib_indicator']
    leveltype = param_lists.translate(settings['level_indicator'])
    level     = settings['level']
    if settings.has_key('tr_ind'): 
        tr_ind = settings['tr_ind']
    #elif source['model'] == Hirlam and 
    else: 
        tr_ind = 0

    print 'In get_grib_plotinfo: getting ',param,settings['level_indicator'],leveltype,level, tr_ind #type(level)

    if grib_index is not None:
        try: # based on all four criteria:
        #if 1:
            selected_grbs=grib_index.select(
                indicatorOfParameter   = param,
                indicatorOfTypeOfLevel = leveltype,
                level                  = level,
                timeRangeIndicator     = tr_ind
                )

        except: # if leveltype is "different"

            print 'Not found: ',param,settings['level_indicator'],leveltype,level, tr_ind 
            return None

        #else:
            selected_grbs=grib_index.select(
                indicatorOfParameter   = param,
                indicatorOfTypeOfLevel = leveltype,
                level                  = level)
            
    else:
        # grib_index may not exist
        print grib_index
        print 'No grib index',param,leveltype,level,tr_ind
        #sys.exit(0)
        return None


    # if there is a selected grib message, retrieve all relevant information
    if selected_grbs:   
        grib_info = {}

        grb = selected_grbs[0]
        grib_info['param']    = param
        grib_info['leveltype'] = leveltype
        grib_info['level']    = grb.level
        if settings.has_key('tr_ind'): grib_info['tr_ind'] = settings['tr_ind']
        else: grib_info['tr_ind'] = 0
        grib_info['values']   = grb.values
        grib_info['date']     = grb.dataDate
        grib_info['time']     = int(grb.dataTime)/100
        grib_info['leadtime'] = max(int(grb.P2),int(grb.startStep) )
        grib_info['latlons']  = grb.latlons()
        if settings['name'] == 'Field': grib_info['name'] = grb.name
        grib_info['model'] = source['name']
        grib_info['source_short'] = source['shortname']

        grib_info = param_lists.set_plotsettings(settings,grib_info)
        print grib_info; sys.exit(0)

        return grib_info
    
        ## no longer executed:
        w         = grb.values
        ddate     = grb.dataDate
        dtime     = grb.dataTime
        leadtime  = int(grb.P2) #int(grb.startStep)
        lats,lons = grb.latlons()

        grib_info = {'param': param,
                     'leveltype': leveltype,
                     'level': grb.level,
                     'name': name,
                     'shortname': shortname,
                     'latlons': (lats,lons),
                     'values': grb.values,
                     'date': grb.dataDate,  # should be starting date of run
                     'time': int(grb.dataTime)/100,  # should be starting time of run
                     'leadtime': max(int(grb.P2),int(grb.startStep) ),
                     'model':source['name'],
                     'source_short':source['shortname'],
                     'plotlevels': levels,
                     'colormap': cmap,
                     'plotcolors':colors,
                     'draw_line': draw_line,
                     'lines': lines
                     }

        return grib_info
                     

    #except: # message not found 
    else:
        print 'No grib message found with ',param,leveltype,level,'time range indicator',tr_ind
        return None

def get_hdf_plotinfo(h5file,settings,source,inlatlon = None):

    import h5tools
    from tools.finddate import datetest

    hdf_info = {}

    print 'In get HDF info'
    
    # background
    hdf_info['model'] = source['model']
    hdf_info['source_short'] = source['shortname']

    h5name = os.path.split(h5file)[1]; print h5name
    datetime,lt = datetest(h5name,pattern = 'ymdh')
    print 'date,time from filename ',h5name, datetime,lt

    # date and time
    hdf_info['date'] = dt.datetime.strftime(datetime,'%Y%m%d')
    hdf_info['time'] = dt.datetime.strftime(datetime,'%H%M')


    #print hdf_info['date'],hdf_info['time']

    # lats,lons
    if inlatlon == None:
        print 'Calculate lats, lons for observation specific projection'
        if 'Radar' in source['model']:
            lats,lons = numpy.array(h5tools.init_h5data(h5file), dtype = float)
            #print 'latlons shape: ', lats.shape, lons.shape
        elif source['model'] == 'MSG':
            lats,lons = numpy.array(h5tools.init_MSG_h5data(h5file), dtype = float)
    else:
        lats,lons = inlatlon

    hdf_info['latlons'] = (lats,lons)

    if 'Radar' in source['model']:
        param = 61
    elif source['model'] == 'MSG':
        param = 71
    else:
        print 'Observations?',source; sys.exit(1)
    hdf_info['param'] = param
    hdf_info['leveltype'] = 'sfc'
    hdf_info['level'] = 0
    hdf_info['leadtime'] = 0
    
    # and the values
    hdf_info['values']    = numpy.array(h5tools.get_h5data(h5file,h5source = source), dtype = float)
    #print 'data shape',values.shape 

    hdf_info = param_lists.set_plotsettings(settings,hdf_info)

    return hdf_info

    # no longer executed

    levels,cmap,colors,draw_line,lines,name,shortname = param_lists.set_plotsettings(settings)
    
    if 1:
        hdf_info  = {'param': param,
                     'grb_indicator': param, # hm, 
                     'leveltype': 'sfc',
                     'level': 0,
                     'name': name,
                     'latlons': (lats,lons),
                     'values': values,
                     'date': date, 
                     'time': time,
                     'leadtime': 0,
                     'model':source['model'],#source['name'],
                     'source_short':source['shortname'],
                     'plotlevels': levels,
                     'colormap': cmap,
                     'plotcolors':colors,
                     'draw_line': draw_line,
                     'lines': lines
                     }

        return hdf_info
