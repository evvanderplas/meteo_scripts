import os,sys,glob,copy,time
import numpy
import sources,param_lists ,settings
import datetime as dt
import tools.finddate
from verif.tools import h5tools
import os, errno

import matplotlib.pyplot
from matplotlib import cm as CM
import gc

sys.path.append('/usr/people/plas/python/tools')
import bmap

from get_plotinfo import *
#print get_plotinfo.__file__; sys.exit(0)

from tools.memory_inspector import report_mem_usage

# create the figure just once, not again and again for each plot
fig  = matplotlib.pyplot.figure()

def suf3(h,suf = 3):
    '''create 0-suffix for integers '''

    if suf == 3:
        if int(h)<10: return '00'+str(h)
        elif int(h)<100: return '0'+str(h)
        elif 1000>int(h)>100:  return '0'+str(int(h)/100)
        elif      int(h)>1000: return str(int(h)/100)
        else: return str(h)
    elif suf == 2:
        if int(h)<10: return '0'+str(h)
        elif 1000>int(h)>100:  return '0'+str(int(h)/100)
        elif      int(h)>1000: return str(int(h)/100)
        else: return str(h)
        

def mkdir_p(path):
    '''Mimic mkdir -p behaviour to create a directory if it does not exist '''

    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise





def set_outdir(plot_info,outdir,domain,datatype='field'):

    if plot_info == None:
        return None

    outdir_dom = os.path.join(outdir,
                              plot_info['source_short'],
                              str(plot_info['date']),
                              suf3(plot_info['time'],suf=2),
                              plot_info['shortname'],
                              str(plot_info['level']),
                              datatype,
                              domain)

    return outdir_dom

def prec_plot(grib_index,grib_index_previous,source,map_list,outdir = './'):

    #print 'oprec_ploty',grib_index

    print 'intprec'
    intprec       = get_grib_plotinfo(grib_index,settings.PCP_i,source)
    print 'accprec',settings.PCP_acc
    accprec       = get_grib_plotinfo(grib_index,settings.PCP_acc,source)
    print 'accprec_prev'
    accprec_prev  = get_grib_plotinfo(grib_index_previous,settings.PCP_acc,source)
    print 'intsnow'
    intsnow       = get_grib_plotinfo(grib_index,settings.SNOW,source,par = 'int')
    print 'accsnow'
    accsnow       = get_grib_plotinfo(grib_index,settings.SNOW,source)
    accsnow_prev  = get_grib_plotinfo(grib_index_previous,settings.SNOW,source)
    intgraup      = get_grib_plotinfo(grib_index,settings.GRAUP,source,par='int')
    accgraup      = get_grib_plotinfo(grib_index,settings.GRAUP,source)
    accgraup_prev = get_grib_plotinfo(grib_index_previous,settings.GRAUP,source)
    #print intprec['values']; sys.exit(0)
    #print accprec['values']; sys.exit(0)
 
    totalprecint = None
    if None not in (intprec,intgraup,intsnow):
        totalprecint  = intprec.copy()
        totalprecint['values']  = 3.6e3 * ( intprec['values'] * 1. + intgraup['values'] * 1. + intsnow['values'] )
    elif intprec is not None:
        totalprecint  = intprec.copy()
        totalprecint['values']  = 3.6e3 * ( intprec['values']) 

    totalaccprec = None
    if accprec:
        totalaccprec  = accprec.copy()
        totalaccsolid = accprec.copy()
        totalaccsolid['name']      = 'Solid precipitation'
        totalaccsolid['shortname'] = 'solid'

        if None not in (accprec_prev,accgraup_prev,accsnow_prev):            
            totalaccprec['values']  = accprec['values'] + accgraup['values'] + accsnow['values']  \
                - accprec_prev['values'] - accgraup_prev['values'] - accsnow_prev['values']
            
            totalaccsolid['values'] = accgraup['values'] + accsnow['values'] \
                - accgraup_prev['values'] - accsnow_prev['values']

        elif accgraup and accsnow:
            totalaccprec['values']  = accprec['values'] + accgraup['values'] + accsnow['values'] 
            totalaccsolid['values'] = accgraup['values'] + accsnow['values'] 
        elif accprec_prev is not None:
            totalaccprec['values']  = accprec['values'] - accprec_prev['values']
            totalaccsolid           = None
        else:
            totalaccprec['values']  = accprec['values']
            totalaccsolid           = None
    else:
        print 'No accprec? ',accprec,grib_index_previous; #sys.exit(0)


    for domain in map_list:
        outdir_dom = set_outdir(intprec,outdir,domain)
        contourplot(intprec,domain = domain,outdir = outdir_dom)

        if 0:
            outdir_dom = set_outdir(accprec,outdir,domain)
            contourplot(accprec,domain = domain,outdir = outdir_dom)
        
        if totalprecint:
            outdir_dom = set_outdir(totalprecint,outdir,domain)
            contourplot(totalprecint,domain = domain,outdir = outdir_dom)

        if totalaccprec:
            outdir_dom = set_outdir(totalaccprec,outdir,domain)
            contourplot(totalaccprec,domain = domain,outdir = outdir_dom)
        
            outdir_dom = set_outdir(totalaccsolid,outdir,domain)
            contourplot(totalaccsolid,domain = domain,outdir = outdir_dom)
        
    return 'Done'

def temp_plot(grib_index,source,map_list,outdir = './'):

    t2m   = get_grib_plotinfo(grib_index,settings.T2M,source)
    try:
        t2m['values'] = t2m['values'] - 273.15
    except:
        print 'No temp output'
        return None

    relh  = get_grib_plotinfo(grib_index,settings.RELH,source)

    for domain in map_list:
        outdir_dom = set_outdir(t2m,outdir,domain)
        contourplot(t2m,domain = domain,outdir = outdir_dom)

    return 0


def cloud_plot(grib_index,source,map_list,outdir = './'):

    clc    = get_grib_plotinfo(grib_index,settings.CLOUD,source)
    brtemp = get_grib_plotinfo(grib_index,settings.BRTEMP,source) # 118 8 sfc 39680 0

    # process:
    try:
        brtemp['values']    = numpy.power( -1 * 17636684 * brtemp['values'], 1/4.) - 273.15
    except:
        print 'No brightness temperature output'

    for domain in map_list:
        outdir_dom = set_outdir(clc,outdir,domain)
        contourplot(clc,domain = domain,outdir = outdir_dom,lsmask_colour = 'red')

        outdir_dom = set_outdir(brtemp,outdir,domain)
        contourplot(brtemp,domain = domain,outdir = outdir_dom,lsmask_colour = 'red')


def vis_plot(grib_index,source,map_list,outdir = './',ftree = True):

    def take_lowest_ml(setting):
        mlset = copy.deepcopy(setting)
        mlset['leveltype'] = 109
        mlset['level']     = 60
        mlset['tr_ind']    = 0

        return mlset
    
    #print settings.PCP_acc
    #print type(settings.PCP_acc)
    #print take_lowest_ml(settings.PCP_acc)
    #print take_lowest_ml(settings.SNOW)
    #sys.exit(1)

    settings.PCP_acc['time'] = 'ml'
    settings.SNOW['time']    = 'ml'
    settings.GRAUP['time']   = 'ml'

    clwat  = get_grib_plotinfo(grib_index,settings.CLWAT,source) 
    tpcp   = get_grib_plotinfo(grib_index,take_lowest_ml(settings.PCP_acc),source,par = 'ml') 
    tsnow  = get_grib_plotinfo(grib_index,take_lowest_ml(settings.SNOW),source,par = 'ml')
    tgraup = get_grib_plotinfo(grib_index,take_lowest_ml(settings.GRAUP),source,par = 'ml')

    if None in (clwat,tpcp,tsnow,tgraup):
        print 'No visibility output ',source['shortname'] 
        return None
    else:
        vis = clwat.copy()

    # # make sure that arguments of power are non-negative, using numpy.select
    dl = 0.
    #dl = 1.e-16
    extinction = (
        144.7 * numpy.power(numpy.select([clwat['values']<0.,clwat['values']  >=0],
                                         [dl,1.2e3*clwat['values']]  ),0.88) + 
        1.1   * numpy.power(numpy.select([tpcp['values']<0.,tpcp['values']>=0],
                                         [dl,1.2e3*tpcp['values']]),0.75) +
        10.4  * numpy.power(numpy.select([tsnow['values']<0.,tsnow['values']>=0],
                                         [dl,1.2e3*tsnow['values']]),0.78) + 
        2.4   * numpy.power(numpy.select([tgraup['values']<0.,tgraup['values']>=0],
                                         [dl,1.2e3*tgraup['values']]),0.78)
        )
    
    vis['values'] = 1000.*3.912 / (extinction + dl) # -ln(0.02) * 1000 / extinction
    
    # get available plot settings
    vis = param_lists.set_plotsettings(settings.VIS,vis)
    
    for domain in map_list:
        if ftree:
            outdir_dom = set_outdir(vis,outdir,domain)
            contourplot(vis, domain = domain,outdir = outdir_dom)
        else:
            contourplot(vis, domain = domain,outdir = outdir)



def windplot(grib_index,source,map_list,outdir = './'):

    #from settings import U10,V10,UGST,VGST

    u10  = get_grib_plotinfo(grib_index,settings.U10,source)
    v10  = get_grib_plotinfo(grib_index,settings.V10,source)
    ugst = get_grib_plotinfo(grib_index,settings.UGST,source)
    vgst = get_grib_plotinfo(grib_index,settings.VGST,source)

    # from m/s to knots: 
    knotfactor = 1852/3600
    knotfactor = 3600./1852.

    totalwind,totalgust = None,None
    try:
        totalwind = u10.copy()
        totalwind = param_lists.set_plotsettings(settings.WIND,totalwind)
    except:
        print 'No wind output'
        return None

    totalwind['values'] = numpy.sqrt( u10['values']**2  + v10['values']**2 ) * knotfactor
    
    try: # only works if gust is in data files
        totalgust = ugst.copy()
        totalgust['values'] = numpy.sqrt( ugst['values']**2 + vgst['values']**2 )
        totalgust = param_lists.set_plotsettings(settings.GUST,totalgust)
    except:
        print 'ugst,vgst (162,163)  not available'

    for domain in map_list:

        if totalwind:
            outdir_dom = set_outdir(totalwind,outdir,domain)
            contourplot(totalwind,domain = domain,outdir = outdir_dom)

            totalwind['shortname'] = 'windvector'
            outdir_dom = set_outdir(totalwind,outdir,domain )
            vectorplot(u10, v10,  domain = domain,outdir = outdir_dom)


        if totalgust:
            outdir_dom = set_outdir(totalgust,outdir,domain)
            contourplot(totalgust,domain = domain,outdir = outdir_dom)



def hdf_plot(source,map_types,latestrun,outdir = './'):

    from verif.tools import h5tools

    if 'Radar' in source['model']:
        hdf_files = sorted( glob.glob( os.path.join(source['dir'],
                                                    source['rexp'])) )

    # as above: choose all (hdf_files) or just the necessary (nhdf_files) files
    for h5file in hdf_files[-18:]:
        print 'opening ',h5file

        # do not put images in separate directory
        if 'Radar' in source['model']:
            path   = os.path.join(outdir,'obs/radar') 
        elif 'MSG' in source['model']:
            path   = os.path.join(outdir,'obs/msg')
        else:
            print 'Do not know where to put the outputfiles of type ',source['model']
            sys.exit(1)

        # check if corresponding file (RAD_NL25_PCP_FC_201211171400.h5, test_RADAR_nl_61_0_2012101903+000.png) 
        # is already there:
        for d in map_types:

            # path to put symbolic link for website
            destdir = os.path.join(outdir,source['shortname'],latestrun[0:8],latestrun[8:10],'pcpi','0','field',d)
            figurefile = os.path.join(path,'obs_'+source['shortname']+'_'+d+'_61_0_'+h5file[-15:-5]+'+000.png')

            if os.path.exists(figurefile):
                print 'plot already there'
                
                # but the symbolic link may not yet be there:
                #destFile = os.path.join(destdir,os.path.split(figurefile)[1])
                #if not os.path.lexists(destFile):
                #    os.symlink(figurefile, destFile)
                    
            else: # plot the contents
                h5info     = get_hdf_plotinfo(h5file,settings.PCP_i,source)
                figurefile = contourplot(h5info,domain = d,
                                         outdir = os.path.join(path,d),
                                         lsmask_colour = 'black')

            # create soft link in findable/viewable dir
            try:
                os.makedirs(destdir)
            except:
                print 'symlink dir already exists: ',destdir
                    
            destFile = os.path.join(destdir,os.path.split(figurefile)[1])
            if not os.path.lexists(destFile):
                os.symlink(figurefile, destFile)
            #sys.exit(1)
                                
                            



def point_extract(grib_index,source,station_dict, ts_object, outdir = './'):

    import tools.interpol as interpol
    import tools.timeseries as ts

    #oper = ts.stationlist()

    for item in settings.point_list:

        param     = item['grib_indicator']
        leveltype = param_lists.translate(item['level_indicator'])
        level     = item['level']

        #if param == 61 and level == 457:
        if param == 11 and leveltype == 'sfc' and level == 2:

            t2m      = get_grib_plotinfo(grib_index,item,source)
            ts_object.addfcdata(t2m)

            if 0: #obsolete?
                for stationid in station_dict.keys():
                    
                    pts   = station_dict[stationid]['points']
                    ll    = station_dict[stationid]['latlon']
                    
                    indic = ((pts[0],pts[2]),(pts[0],pts[3]),(pts[1],pts[3]),(pts[1],pts[2]))
                    value = interpol.interp(t2m['values'],t2m['latlons'],ll,indic)

                    #station_dict[stationid]['value'] = interpol.interp(t2m['values'],t2m['latlons'],ll,indic)
                    if 1:
                        #print ll,pts
                        #print station_dict[stationid]
                        for s in indic: print 'at corners: ', t2m['values'][s]
                    
                        print station_dict[stationid]['value']
                    #pass
        


def contourplot(grib_info,domain = 'nl',outdir = './',
                lsmask_colour = 'black',exarg = 'both',
                v=0):
    
    global fig

    if grib_info == None:
        print 'No plot info'
        return None


    #fig  = matplotlib.pyplot.figure()
    ax   = fig.add_subplot(1,1,1)

    t1 = time.clock()
    if v: print 'contourplot: dressing up',time.clock() - t1
    plot_map = bmap.myMap(domain = domain, modelname = grib_info['model'])
    plot_map.xymap(latlons = grib_info['latlons'], domain = domain, modelname = grib_info['model'])
    plot_map.set_axes(ax)
    plot_map.dress_up_map(domain = domain,lsmask_colour = lsmask_colour)
    if v: print 'done dressing up',time.clock() - t1
    
    diy = True # switch to use private plotting (True) or testmap in myMap instance (False)
    if diy:
        if grib_info['plotlevels'] == None:
            gmax,gmin = grib_info['values'].max(),grib_info['values'].min()
            grib_info['plotlevels'] = numpy.arange(gmin,gmax,(gmax-gmin)/10.)
            print gmin,gmax,grib_info['plotlevels']
        if grib_info['colormap'] == None and grib_info['plotcolors'] == None: 
            cmap = CM.get_cmap('jet')
        elif grib_info['colormap'] is not None:
            cmap = grib_info['colormap']
            cmap.set_over(color='k',alpha = None)
            grib_info['plotcolors'] = None

            # for gray-scale pictures, do not extend the colorbar beyond [0,1]
            #for item in grib_info.keys(): print item, grib_info[item]
            #if grib_info['grib_indicator'] == 71: # does not always work?
            if grib_info['param'] == 71: 
                #if cmap in (CM.get_cmap('gray'),CM.get_cmap('gray_r'),CM.get_cmap('binary_r')):
                exarg = 'neither'

        elif grib_info['plotcolors'] is not None:
            cmap = None
        else:
            pass

        # validtime:
        # for observations, make sure leadtime is 0 (were inserted into arbitrary grib message)
        if grib_info['model'] in ('Radar NL','Radar EUR','MSG'):
            grib_info['leadtime'] = 0 

        if int(grib_info['time']) > 100: # could be refined to insert minutes for example
            grib_info['time'] = int(grib_info['time'])/100

        validtime = dt.datetime.strptime(str(grib_info['date'])[0:8],'%Y%m%d') + \
            dt.timedelta(hours = int(grib_info['time']) ) + \
            dt.timedelta(hours = int(grib_info['leadtime']) )
        validdtg = validtime.strftime('%Y-%m-%d, %H')

        #print str(grib_info['date']),str(grib_info['time']),dt.datetime.strptime(str(grib_info['date'])[0:8],'%Y%m%d'),validtime, validdtg; #sys.exit(1)

        #title
        ax.set_title( grib_info['model']+': '+grib_info['name']+' \n(parameter '+str(grib_info['param'])+\
                          ', level '+str(grib_info['level'])+\
                          ') \n at '+str(grib_info['date'])+suf3(grib_info['time'],suf=2)+\
                          '+'+suf3(grib_info['leadtime'])+': '+str(validtime)+' UTC') 


        # actual plots
        if v: print 'actual plot',time.clock() - t1
        if v: print 'x,y,data shape',plot_map.x.shape,plot_map.y.shape,grib_info['values'].shape
        pcont = plot_map.bmap.contourf(plot_map.x,plot_map.y,grib_info['values'],
                                       grib_info['plotlevels'],
                                       cmap=cmap,colors = grib_info['plotcolors'],extend=exarg,ax=ax)

        if grib_info['draw_line']: 
            if v: print 'drawing contour lines: ',
            if grib_info.has_key('lines'): 
                lines = grib_info['lines']
                print lines
            else:
                print 'contourplot draw_line: is this an option?'
                lines = grib_info['plotlevels'] # probably redundant
            plcont = plot_map.bmap.contour(plot_map.x,plot_map.y,\
                                               grib_info['values'],lines,colors='black',ax=ax)

        # colorbar:
        cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.2f') #format='%.1e') 
        cb.set_label(grib_info['name'],fontsize=9)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(9)
        if grib_info['draw_line']: 
            cb.add_lines(plcont)

        # naming and saving
        # include source,domain, and parameter/level, and starting date + leadtime
        if v: print 'getting plot title',time.clock() - t1        
        mkdir_p(outdir)

        plotfile = os.path.join(outdir,
                                'test_'+grib_info['source_short']+'_'+str(domain)+'_'\
                                +str(grib_info['param'])+'_'+str(grib_info['level'])+'_'\
                                +str(grib_info['date'])+suf3(grib_info['time'],suf=2)+'+'+\
                                suf3(grib_info['leadtime'])+'.png')
        if v: print 'saving plot',time.clock() - t1        
        fig.savefig(plotfile, dpi = 100) #dpi = 300)

        cmd = 'convert -trim '+plotfile+' '+plotfile
        #os.system(cmd)
        print 'created plot: '+plotfile,time.clock() - t1


        textfile = os.path.join(outdir,'README.txt')
        if os.path.isfile(textfile):
            pass
        else:
            f = open(textfile,'wt')
            f.write('Plot generated by fabriek.py (plottypes.contourplot()). \n\nContact plas@knmi.nl')
            f.close()

        #del(cbcoll)
        del(ax)

        fig.clf()
        gc.collect() # force garbage collection

        return plotfile

        
def vectorplot(grib_info_x, grib_info_y, domain = 'nl', outdir = './', grib_info = {}):
                
    '''A quiver plot based on two fields, perhaps add an optional underlying field '''

    global fig
    
    # fig  = matplotlib.pyplot.figure()
    ax   = fig.add_subplot(1,1,1)

    plot_map = bmap.myMap(domain = domain, modelname = grib_info_x['model'])
    plot_map.xymap(latlons = grib_info_x['latlons'], domain = domain, modelname = grib_info_x['model'])
    plot_map.set_axes(ax)
    plot_map.dress_up_map(domain = domain)
     
    if grib_info.has_key('parameter'):
        print 'huh?',grib_info
        if grib_info['colormap'] == None and  grib_info['plotcolors'] == None: 
            cmap = CM.get_cmap('jet')
        elif grib_info['colormap'] is not None:
            cmap = grib_info['colormap']
            cmap.set_over(color='k',alpha = None)
            grib_info['plotcolors'] = None
        else:
            pass
        
    #title
    ax.set_title( grib_info_x['model']+': '+grib_info_x['name']+' (parameter '+str(grib_info_x['param'])+\
                      ', level '+str(grib_info_x['level'])+\
                      ') \n at '+str(grib_info_x['date'])+suf3(grib_info_x['time'],suf=2)+'+'+suf3(grib_info_x['leadtime']) 
                  )

    #thinning factor:
    if domain == 'nl' and 'Harmonie' in grib_info_x['model']:
        thf = 8
    elif domain == 'eur':
        thf = 20
    elif domain in ('rot','ijs') and 'Harmonie' in grib_info_x['model'] :
        thf = 1
    else:
        thf = 8
        

    # actual plots
    arrows = plot_map.bmap.quiver(plot_map.x[::thf,::thf],plot_map.y[::thf,::thf],
                                  grib_info_x['values'][::thf,::thf],grib_info_y['values'][::thf,::thf],
                                  scale=400,
                                  color = 'red',
                                  ax=ax)

 
    if 0: # maybe to plot second
        # colorbar:
        cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.1e') 
        cb.set_label(grib_info['name'],fontsize=9)

    # naming and saving
    mkdir_p(outdir)

    plotfile = os.path.join(outdir,
                            'test_vec_'+grib_info_x['source_short']+'_'+str(domain)+'_'\
                                +str(grib_info_x['param'])+'_'+str(grib_info_x['level'])+'_'\
                                +str(grib_info_x['date'])+suf3(grib_info_x['time'],suf=2)+'+'+suf3(grib_info_x['leadtime'])+'.png')

    sf = fig.savefig(plotfile)
    print 'created plot: '+plotfile

    textfile = os.path.join(outdir,'README.txt')
    if os.path.isfile(textfile):
        pass
    else:
        f = open(textfile,'wt')
        f.write('Plot generated by fabriek.py (plottypes.vectorplot()). \n\nContact plas@knmi.nl')
        f.close()

        
    #print dir(arrows)
    #sys.exit(1)
    #clean up
    #for coll in arrows.collections:
    #    coll.remove()
    
    try:
        for coll in pcont.collections:
            coll.remove()

        for coll in plcont.collections:
            coll.remove()

        # delete colorbar
        cbcoll = fig.get_axes()[1]
        cbcoll.collections = []
        fig.delaxes(cbcoll)
        fig.subplots_adjust(right=0.90)

    except:
        pass


