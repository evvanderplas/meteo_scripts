import os,sys,glob,copy,time
import numpy
import sources,param_lists ,settings
import datetime as dt
#import tools.finddate
#from verif.tools import h5tools
import os, errno
import matplotlib.pyplot
import gc

sys.path.append('/usr/people/plas/python/tools')
import bmap

from get_plotinfo import *

from tools.memory_inspector import report_mem_usage

# create the figure just once, not again and again for each plot
fig  = matplotlib.pyplot.figure()

def mkdir_p(path):
    '''Mimic mkdir -p behaviour to create a directory if it does not exist '''

    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def get_inf(gfile,setting,source):
    
    '''
    Use python/tools plotinfo to extract information from grib to plot it, manipulate it

    @ 20140129: 
    def read_from_grib(self,gribfile,par = 61,level = 0,leveltype='sfc',TR=0,name=None,case='test',model = 'Harmonie')
    def read_from_h5(self,h5file,datapath = 'image1/image_data',metapath='overview',name='precipitation')
    '''

    import plotinfo

    if gfile is None: return None

    try: # is grib?:
        if setting.has_key('leveltype'):
            leveltype = setting['leveltype']
        elif setting.has_key('level_indicator'):
            leveltype = param_lists.translate(setting['level_indicator'])
        else:
            leveltype = 'sfc'
            
        if not setting.has_key('tr_ind'):
            tr = 0
        else:
            tr = setting['tr_ind']
    
        # check which setting to use: hir for hirlam, new As of 37h1.2 or default
        if setting.has_key('grib_indicator_hir') and source['model'] == 'Hirlam':
            par   = setting['grib_indicator_hir']
            level = setting['level_hir']
            tr    = setting['tr_ind_hir']
        elif setting.has_key('grib_indicator_n') and int(source['version']) > 3711:
            par   = setting['grib_indicator_n']
            level = setting['level_n']
            tr    = setting['tr_ind_n']
        else:
            par   = setting['grib_indicator']
            level = setting['level']
            #tr    = setting['tr_ind']
    
    except KeyError:
        print 'plot hdf file:'
        pass
    
    inf = plotinfo.plotinfo(model='init')

    print 'Test hdf5 ',gfile
    if type(gfile) == str:
        if os.path.splitext(gfile)[1] in ('.h5','.hdf5'):

            ## EvdP this is not elegant: either sources should have a location where to find time (so metapath,dtattr in sources.py), or 
            ## a try/except loop should be included to somehow circumvent this issue
            if source.has_key('dtattr'): 
                inf.read_from_h5(gfile,datapath=source['datapath'],metapath=source['metapath'],dtattr = source['dtattr'],model=source['shortname'],name=source['name'])
            else:
                inf.read_from_h5(gfile,datapath=source['datapath'],metapath=source['metapath'],model=source['shortname'],name=source['name'])
        else:
            inf.read_from_grib(gfile,par=par,leveltype=leveltype,level=level,TR=tr,model=source['name'])
    else:
        inf.read_from_grib(gfile,par=par,leveltype=leveltype,level=level,TR=tr,model=source['name'])
    #inf.get_settings(preset = setting['shortname'])

    # incorporate name conventions, colors etc from settings
    for k,v in setting.iteritems():
        #print k,v
        inf.__dict__[k] = v  
    inf.source = source['shortname']
    
    return inf


def set_outdir(plot_info,outdir,domain,datatype='field'):

    if plot_info == None:
        return None



    if type(plot_info) == dict:
        srcname         = plot_info['source_short']
        pdate,ptime     = str(plot_info['date']),str(plot_info['time']).zfill(2)
        shortname,level = plot_info['shortname'],str(plot_info['level'])
    else:
        srcname         = plot_info.source
        pdate,ptime     = str(plot_info.date),str(plot_info.time).zfill(2)
        shortname,level = plot_info.shortname,str(plot_info.level)

    #print 'outdir',outdir,srcname,pdate,ptime,shortname,level,datatype,domain
    outdir_dom = os.path.join(outdir,
                              srcname,
                              pdate,ptime,
                              shortname,level,
                              datatype,
                              domain)

    # Create directory if it does not exist:
    mkdir_p(outdir_dom)

    return outdir_dom

def prec_plot(grib_index,grib_index_previous,source,map_list,outdir = './'):

    #print 'oprec_ploty',grib_index
        
    accprec       = get_inf(grib_index,settings.PCP_acc,source)
    intprec       = get_inf(grib_index,settings.PCP_i,source)
    accprec_prev  = get_inf(grib_index_previous,settings.PCP_acc,source)

    intsnow       = get_inf(grib_index,settings.SNOW_i,source)
    accsnow       = get_inf(grib_index,settings.SNOW,source)
    accsnow_prev  = get_inf(grib_index_previous,settings.SNOW,source)

    intgraup      = get_inf(grib_index,settings.GRAUP_i,source)
    accgraup      = get_inf(grib_index,settings.GRAUP,source)
    accgraup_prev = get_inf(grib_index_previous,settings.GRAUP,source)
    
    #print intprec.values.shape, intsnow.values.shape, intgraup.values.shape
    totalprecint = 3600. * (intprec + intsnow + intgraup)
    
    totalaccprec  = accprec - accprec_prev + accsnow - accsnow_prev + accgraup - accgraup_prev
    totalaccsolid = accsnow - accsnow_prev + accgraup - accgraup_prev
    totalaccsolid.shortname = 'solid' 

    for domain in map_list:
        outdir_dom = set_outdir(intprec,outdir,domain)
        #print outdir_dom; sys.exit(0)
        totalprecint.plot(domain=domain,out = outdir_dom)

        outdir_dom = set_outdir(accprec,outdir,domain)
        totalaccprec.plot(domain=domain,out = outdir_dom)

        outdir_dom = set_outdir(totalaccsolid,outdir,domain)
        totalaccsolid.plot(domain=domain,out = outdir_dom)
        
    return 'Done'

def temp_plot(grib_index,source,map_list,outdir = './'):

    t2m   = get_inf(grib_index,settings.T2M,source)

    try:
        if t2m.values.mean() > 200.:
            t2m = t2m - 273.15
    except:
        print 'No temp output'
        return None

    relh  = get_inf(grib_index,settings.RELH,source)
    q2m   = get_inf(grib_index,settings.Q2M,source)

    for domain in map_list:
        outdir_dom = set_outdir(t2m,outdir,domain)
        t2m.plot(domain = domain,out = outdir_dom)
        outdir_dom = set_outdir(q2m,outdir,domain)
        q2m.plot(domain = domain,out = outdir_dom)
        outdir_dom = set_outdir(relh,outdir,domain)
        relh.plot(domain = domain,out = outdir_dom)

    return 'Done'


def cloud_plot(grib_index,source,map_list,outdir = './'):

    clc    = get_inf(grib_index,settings.CLOUD,source)
    brtemp = get_inf(grib_index,settings.BRTEMP,source) # 118 8 sfc 39680 0

    #return brtemp
    # process:
    try:
        brtemp = (-1 * 17636684 * brtemp)**( 1/4.) - 273.15
    except:
        print 'No brightness temperature output'

    for domain in map_list:
        outdir_dom = set_outdir(clc,outdir,domain)
        clc.plot(domain = domain,out = outdir_dom,lsmask_colour = 'red')

        outdir_dom = set_outdir(brtemp,outdir,domain)
        brtemp.plot(domain = domain,out = outdir_dom,lsmask_colour = 'red')


def vis_plot(grib_index,source,map_list,outdir = './',ftree = True):

    def take_lowest_ml(setting):

        # do not change setting itself:
        mlset = copy.deepcopy(setting)
        
        ## on model levels different grib indicator, sigh...
        if   mlset['grib_indicator'] == 61: mlset['grib_indicator'] = 62 # rain, on ml != sfc...
        elif mlset['grib_indicator'] == 62: mlset['grib_indicator'] = 79
        elif mlset['grib_indicator'] == 63: mlset['grib_indicator'] = 201

        ## use lowest model level:
        mlset['leveltype'] = 109 #'ml' #109
        mlset['level']     = 60
        mlset['tr_ind']    = 0

        return mlset

    clwat  = get_inf(grib_index,settings.CLWAT,source) 
    tpcp   = get_inf(grib_index,take_lowest_ml(settings.PCP_acc),source) 
    tsnow  = get_inf(grib_index,take_lowest_ml(settings.SNOW),source)
    tgraup = get_inf(grib_index,take_lowest_ml(settings.GRAUP),source)

    for hydrom in (clwat,tpcp,tsnow,tgraup):
        if hydrom.values == []:
            print 'No visibility output ',source['shortname'] 
            return None
    else:
        vis = clwat.copy()

    # # make sure that arguments of power are non-negative, using numpy.select
    dl = 0.
    #dl = 1.e-10
    extinction = (
        144.7 * numpy.power(numpy.select([clwat.values  < 0.,clwat.values  >= 0],
                                         [dl,1.2e3*clwat.values]  ),0.88) + 
        1.1   * numpy.power(numpy.select([tpcp.values   < 0.,tpcp.values   >= 0],
                                         [dl,1.2e3*tpcp.values]),0.75) +
        10.4  * numpy.power(numpy.select([tsnow.values  < 0.,tsnow.values  >= 0],
                                         [dl,1.2e3*tsnow.values]),0.78) + 
        2.4   * numpy.power(numpy.select([tgraup.values < 0.,tgraup.values >= 0],
                                         [dl,1.2e3*tgraup.values]),0.78)
        )
    
    vis.values = 1000.*3.912 / (extinction + dl) # -ln(0.02) * 1000 / extinction
    
    # get available plot settings
    for k,v in settings.VIS.iteritems(): 
        vis.__dict__[k] = v
        #print k,v

    # plot
    for domain in map_list:
        outdir_dom = set_outdir(vis,outdir,domain)
        vis.plot(domain = domain,out = outdir_dom)


def windplot(grib_index,source,map_list,outdir = './'):

    #from settings import U10,V10,UGST,VGST

    u10  = get_inf(grib_index,settings.U10,source)
    v10  = get_inf(grib_index,settings.V10,source)
    ugst = get_inf(grib_index,settings.UGST,source)
    vgst = get_inf(grib_index,settings.VGST,source)

    # from m/s to knots: 
    knotfactor = 1852/3600
    knotfactor = 3600./1852.

    totalwind = ( u10*u10  + v10*v10 )**(0.5) * knotfactor
    for k,v in settings.WIND.iteritems(): totalwind.__dict__[k] = v

    totalgust = ( ugst*ugst  + vgst*vgst)**(0.5) 
    for k,v in settings.GUST.iteritems(): totalgust.__dict__[k] = v
    
    for domain in map_list:
        outdir_dom = set_outdir(totalwind,outdir,domain)
        totalwind.plot(domain = domain,out = outdir_dom)

        outdir_dom = set_outdir(totalgust,outdir,domain)
        totalgust.plot(domain = domain,out = outdir_dom)



def hdf_plot(source,map_types,latestrun,outdir = './'):

    #from verif.tools import h5tools

    try:
        ymd = latestrun.strftime('%Y%m%d')
        st  = latestrun.strftime('%H')
    except:
        try:
            
            ymd = str(latestrun)[0:8]
            st  = str(latestrun)[8:10]
            latestrun = dt.datetime.strptime(latestrun[0:10],'%Y%m%d%H')
        except:
            print 'date not recognized',latestrun
    
    #sys.exit(0)

    if 'Radar' in source['model']:
        hdf_files = sorted( glob.glob( os.path.join(source['dir'],
                                                    source['rexp'])) )

    # as above: choose all (hdf_files) or just the necessary (nhdf_files) files
    for dl in range(24,0,-1):
        h = latestrun - dt.timedelta(hours=dl)
        if 'msg' in source['name'].lower() and 'ir' in source['name'].lower():
            setting = settings.BRTEMP
            h5file = os.path.join(source['dir'],source['fformat'].format(y = h.strftime('%Y'),m = h.strftime('%m'),d = h.strftime('%d'),h = h.strftime('%H')))
            outpath   = os.path.join(outdir,'obs/msg')
        elif 'msg' in source['name'].lower() and 'cm' in source['name'].lower():
            setting = settings.CLOUD
            h5file = os.path.join(source['dir'].format(y = h.strftime('%Y'),m = h.strftime('%m'),d = h.strftime('%d')),source['fformat'].format(dtgh = h.strftime('%Y%m%d%H')))
            outpath   = os.path.join(outdir,'obs/msg')
        else:
            setting = settings.PCP_i
            h5file = os.path.join(source['dir'],source['fformat'].format(dtgh = h.strftime('%Y%m%d%H')))
            outpath   = os.path.join(outdir,'obs/radar') 
        if not os.path.exists(h5file): break

        #for h5file in hdf_files[-18:]:
        print 'opening ',h5file
        #for k,v in setting.iteritems(): print k,v


        # check if corresponding file (RAD_NL25_PCP_FC_201211171400.h5, test_RADAR_nl_61_0_2012101903+000.png) 
        # is already there:
        for d in map_types:

            # path to put symbolic link for website
            destdir = os.path.join(outdir,source['shortname'],ymd,st,setting['shortname'],'0','field',d)
            # all obs fields into one column:
            destdir = os.path.join(outdir,'obs',ymd,st,setting['shortname'],'0','field',d)
            destfile = 'obs_{name}_{dmap}_0_{ymd}{st}+000.png'.format(name=source['shortname'],dmap=d,ymd=h.strftime('%Y%m%d'),st=h.strftime('%H'))
            figurefile = os.path.join(outpath,destfile)

            if os.path.exists(figurefile):
                print 'plot already there',figurefile
                
                # but the symbolic link may not yet be there:
                #destFile = os.path.join(destdir,os.path.split(figurefile)[1])
                #if not os.path.lexists(destFile):
                #    os.symlink(figurefile, destFile)
                    
            else: # plot the contents
                sgribsetting = settings.PCP_acc
                if d == 'eur':
                    sfile = '/nobackup_1/users/plas/python/fabriek/fabriek/samplegrib_LA.grb' 
                    sample = get_inf(sfile,sgribsetting,sources.BULL_OPER)
                elif d == 'nl':
                    sfile = '/nobackup_1/users/plas/python/fabriek/fabriek/samplegrib_SA.grb' 
                    sample = get_inf(sfile,sgribsetting,sources.BULL_OPER_SA)            
                else:
                    print 'Not sure, resample to european area'
                    sample = get_inf('samplegrib_LA',sgribsetting,sources.BULL_OPER)

                h5info     = get_inf(h5file,setting,source); print h5info.date
                if h5info.latlons: h5info.resample(sample,nodata = 65535);      print h5info.date
                figurefile = h5info.plot(domain = d,out = figurefile,lsmask_colour = 'red') #'black')
                if figurefile == None: 
                    print 'plotting failed'
                    return None
                

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
                          ') \n at '+str(grib_info['date'])+str(grib_info['time']).zfill(2)+\
                          '+'+str(grib_info['leadtime']).zfill(3)+': '+str(validtime)+' UTC') 


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
                                +str(grib_info['date'])+str(grib_info['time']).zfill(2)+'+'+\
                                str(grib_info['leadtime']).zfill(3)+'.png')
        if v: print 'saving plot',time.clock() - t1        
        fig.savefig(plotfile, dpi = 100) #dpi = 300)

        cmd = 'convert -trim '+plotfile+' '+plotfile
        os.system(cmd)
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
                      ') \n at '+str(grib_info_x['date'])+str(grib_info_x['time']).zfill(2)+'+'+str(grib_info_x['leadtime']).zfill(3) 
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
                                +str(grib_info_x['date'])+str(grib_info_x['time']).zfill(2)+'+'+str(grib_info_x['leadtime']).zfill(3)+'.png')

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


