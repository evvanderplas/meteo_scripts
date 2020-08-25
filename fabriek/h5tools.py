#!/usr/bin/env python

import h5py
import os,sys,glob,numpy
import datetime as dt

def generate_radar_latlons(dim0,dim1,ll_lon,ll_lat,ur_lon,ur_lat):

    '''Generate lats,lons for polar stereographically projected radar from map characteristics '''
    
    # import basemap:
    from mpl_toolkits.basemap import Basemap # map functionality


    maph = Basemap(projection = 'stere',
                   lat_0 = 90.,  lon_0 = 0.,
                   lat_ts = 60., 
                   llcrnrlon = ll_lon, llcrnrlat = ll_lat,
                   urcrnrlon = ur_lon, urcrnrlat = ur_lat 
                   )

    lons,lats = maph.makegrid(dim1,dim0,returnxy=False)

    return lats,lons

def generate_MSG_latlons(dim0,dim1,test = False):

    '''Generate lats,lons for MSG geostationary satellite image from map characteristics '''
    
    # import basemap:
    from mpl_toolkits.basemap import Basemap # map functionality


    if test:
        # test: first make a projection of whole earth
        m1 = Basemap(projection='geos',lon_0 = 0.0, resolution = None)
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')

        print 'llx,y ',m1.llcrnrx,m1.llcrnry # 0 0
        print 'urx,y ',m1.urcrnrx,m1.urcrnry # 10857955.0487 10857955.0487

    # from hdf5 file:
    #proj string: +proj=geos +a=6378169.0 +b=6356583.8 +lon_0=0.0 +h=35785831.0
    x_lr =  5565748.2275008615
    x_ul = -5568748.630858005
    y_lr =  2655356.9710718393
    y_ul =  5436730.883143699

    maph = Basemap(projection = 'geos',
                   rsphere=(6378169.00,6356583.8),
                   #resolution='l', #None,
                   #area_thresh=10000.,
                   lon_0=0.0,
                   llcrnrx = x_ul,llcrnry = y_lr,
                   urcrnrx = x_lr,urcrnry = y_ul
                   #llcrnrx= -m1.urcrnrx/2.,
                   #llcrnry= m1.llcrnry, #0.0,
                   #urcrnrx= m1.urcrnrx/2.,
                   #urcrnry= m1.urcrnry/2.
                   ,satellite_height=35785831.0
                   )

    if test:
        xa,ya   = maph(4.878867,52.375345) # location of NL148
        xc,yc   = maph(0.,45.)             # location sw in france
        print xc,yc
        maph.drawcoastlines(color = 'w')
        maph.drawcountries(color = 'w')
        maph.fillcontinents(color='coral',lake_color='aqua')
        maph.drawmapboundary(fill_color='aqua')
        maph.plot([xa,xc],[ya,yc],'o',color='white')
        plt.show() # picture looks exactly like hdfview image!

    xtot = x_lr - x_ul
    ytot = y_ul - y_lr
    xs = numpy.arange(x_ul + xtot/(2*dim1),x_lr,xtot/dim1)
    ys = numpy.arange(y_lr + ytot/(2*dim0),y_ul,ytot/dim0)
    print len(xs),len(ys),min(xs),x_ul,min(ys),max(ys) #; sys.exit(1)

    x  = numpy.array([xs for i in range(dim0)])
    y  = numpy.array([ys for i in range(dim1)]).transpose()

    print x.shape, y.shape

    #lonpt, latpt  = maph(x,y,inverse=True); #print lonpt,latpt
    lons,lats,x,y = maph.makegrid(dim1,dim1,returnxy=True) #False
    #lons,lats,x,y = maph.makegrid(dim0,dim1,returnxy=True) #False
    print '4kant: ',lons.shape, lons[100,1856],lons[-100,1856], lons[1856,100],lons[1856,-100]
    print '4kant: ',lats.shape, lats[100,1856],lats[-100,1856], lats[1856,100],lats[1856,-100]#; sys.exit(1)


    #print lonpt.shape, latpt.shape
    a,b = 1856,1856
    print 'zentrum       ', lons[a,b], lats[a,b],x[a,b],y[a,b]
    #print 'zentrum selbst', lonpt[464,1856],latpt[464,1856],xs[464,1856],ys[464,1856]
    if test:
        #lons,lats,x,y = maph.makegrid(dim0,dim1,returnxy=True) #False
        #print 'eh', min(lons[:,10]),max(lons[:,10]),min(lats[1800,:]), max(lats[1800,:])
        #print 'eh', min(x[10,:]),max(x[10,:]),min(y[:,1800]), max(y[1800,:]); #sys.exit(1)

        print 'eh', min(lons[10,:]),max(lons[10,:]),min(lats[:,1800]), max(lats[:,1800])
        print 'eh', min(x[10,:]),max(x[10,:]),min(y[:,1800]), max(y[:,1800]); sys.exit(1)
    #lons,lats = m1.makegrid(dim1,dim0,returnxy=False) # makes no difference!

    #return latpt,lonpt
    return lats[-928:,:],lons[-928:,:]

def init_h5data(file):
    
    '''Initialise map for hdf5 radar data '''

    import h5py

    f=h5py.File(file,'r')

    # once you know where the data is, fill in 'path':
    path = '/image1/image_data'
    data = f[path]
    w = data[:,:] #/100
    (dim0,dim1) = w.shape
    f.close()

    if dim0 == 765:

        (ll_lat,ll_lon) = (49.362,  0. )
        (ur_lat,ur_lon) = (55.389, 10.856)
                    
    elif dim0 == 256:

        (ll_lat,ll_lon) = (49.769,  0.)
        (ur_lat,ur_lon) = (54.818,  9.743)
        
    elif dim0 == 512:

        (ll_lat,ll_lon) = (41.937, -9.271)
        (ur_lat,ur_lon) = (58.089, 20.453)
            
    #print 'dimensions, ll corners',dim0,dim1,ll_lon,ll_lat,ur_lon,ur_lat
    lats,lons = generate_radar_latlons(dim0,dim1,ll_lon,ll_lat,ur_lon,ur_lat)

    return lats,lons

def init_MSG_h5data(file):
    
    '''Initialise map for hdf5 MSG data '''

    import h5py

    f=h5py.File(file,'r')

    # once you know where the data is, fill in 'path':
    path = '/CMa'
    data = f[path]
    w = data[:,:] #/100
    (dim0,dim1) = w.shape
    f.close()

    print dim0,dim1
    #sys.exit(1)
            
    #print 'dimensions, ll corners',dim0,dim1,ll_lon,ll_lat,ur_lon,ur_lat
    lats,lons = generate_MSG_latlons(dim0,dim1) #,test = True)

    print w[450,1856]
    #print lats.shape,lats[450,1856],lats[920,1856],lats[9,1856],lats[9,1650],lats[9,2250]
    #print lons.shape,lons[450,1856],lons[920,1856],lons[9,1856],lons[9,1650],lons[9,2050]


    #print lons[10,1450:1750]
    #print lats[10:500,1750]
    #sys.exit(1)
    return lats,lons

def get_h5data(hfile,h5source = None):
    '''
    Extract hdf5 radar data from standard path /image1/image_data 
    
    Returns rain rate if the name of the source is "Radar" or None,
    returns cloud cover (0, 0.1 or 1) from cloud mask SAF for "MSG"
    '''
    
    import h5py
    import numpy

    try:
        f     = h5py.File(hfile,'r')
    except IOError:
        print 'File cannot be opened: ',hfile
        return None

    fname = os.path.split(hfile)[1]


    if h5source == None or 'Radar' in h5source['model']: # default to radar

        # once you know where the data is, fill in 'path':
        path = '/image1/image_data'

        data = f[path]
        w = data[:,:] 
        f.close()

        nodata = 65535
        if '_RA' in fname:
            print 'accumulated rainrate: '#,max(w)
            nodata_crit = [w[::-1,:]<nodata,w[::-1,:] == nodata]
            #nodata_do   = [0.01 * w[::-1,:],w[::-1,:]]
            nodata_do   = [0.01 * w[::-1,:],nodata] #0]
        else:
            print 'radar reflectivity: '#,max(w)
            nodata_crit = [w[::-1,:]<255,w[::-1,:] == 255]
            nodata_do   = [numpy.power(10,(w[::-1,:] - 109.)/32.),nodata]

        wn = numpy.select(nodata_crit,nodata_do)

    elif h5source['model'] == 'MSG':

        # once you know where the data is, fill in 'path':
        path = '/CMa'

        data = f[path]
        w   = data[:,:]  
        # wn.where 1 = cloud free, 2 is partly cloudy ("contaminated"), 3 = fully clouded, 4 snow/ice, 5 unprocessed 
        # http://www.nwcsaf.org/HD/MainNS.jsp
        nodata = 65535
        wn = numpy.select([w<1,w == 1,w==2,w==3,w>3],[nodata,0,0.1,1.,nodata])
        f.close()

        wn   = wn[::-1,:] 

    #print 'now get: ',wn.max(),wn.mean(),wn.min()
    #print 'hist:',numpy.histogram(wn,bins=[0,1,2,3,4,5])
    return wn

def resample_radar(raddata,radlats,radlons,newlats,newlons,nodata = 65535):

    import pyresample

    polar_def  = pyresample.geometry.GridDefinition(lons=radlons, lats=radlats)
    new_def    = pyresample.geometry.GridDefinition(lons=newlons, lats=newlats)


    msg_con_nn     = pyresample.image.ImageContainerNearest(raddata, polar_def, 
                                                            radius_of_influence=20000, fill_value = nodata)
    area_con_nn    = msg_con_nn.resample(new_def)
    result_data_nn = area_con_nn.image_data

    return result_data_nn

def accumulate_radar(h5files,gribfiles,outdir = './',accstep=1,acc = True, start = 0):

    '''Finds reflectivity(!) data in a hdf5 file, 
    resamples it to fc grid and puts it in a grib message  '''

    if h5files == [] or gribfiles == []:
        print 'Accumulate radar: Files not found'
        return 0

    import find_date,gribtools,config

    config.mkdir_p(outdir)

    glats,glons = gribtools.get_latlon(gribfiles[0])
    rlats,rlons = init_h5data(h5files[0])

    def tofile(gfile,date,data,step,param=61,level=457):

        resampled_data = resample_radar(data,rlats,rlons,glats,glons)

        #print 'now resampled: ',resampled_data.max(), resampled_data.mean(), resampled_data.min()
        dtg     = dt.datetime.strftime(date,'%Y%m%d%H')
        newfile = os.path.join(outdir,'rad_'+str(int(step))+'_'+str(dtg)+'+000.grib')
        gribtools.insert_grib(gfile, newfile, resampled_data,parameter = param, level = level)
        #gribtools.insert_grib(gfile, newfile, resampled_data,parameter = 0, level = level)

        return newfile

    def accumulate(gfile,h5files):

        w = None #? 

        for afile in h5files:

            hdate = find_date.datetest(afile,pattern='hdf')[0]
            print 'files:',os.path.split(gfile)[1],os.path.split(afile)[1]
            print gdate, hdate, acctime.seconds/3.6e3,gdate - acctime < hdate <= gdate
            
            if gdate - acctime < hdate <= gdate: # obs time between fctime and an accumulation time before

                wa    = get_h5data(afile) 
                if acc:
                    try:
                        w = w + wa
                        nr += 1
                    except: # first
                        w  = wa
                        nr = 1
                    print 'now: nr, max, mean, min: ',nr, w.max(), w.mean(), w.min()
                else: # for the RACfiles: 3h accumulation every hour available!
                    #print 'now tot: ',wa.max(), wa.min()
                    if gdate - dt.timedelta(hours=1) < hdate <= gdate:
                        accfile = tofile(gfile,gdate,wa,acctime.seconds/3.6e3)
                        #radfiles.append( accfile )
                        return accfile
                    else:
                        pass

            elif hdate > gdate or hdate <= gdate - acctime : # time to write to file
                if acc:
                    #try:
                    if w not in [[],None]: # check if data to write to file is not None
                        accfile = tofile(gfile,gdate,w/nr,acctime.seconds/3.6e3)

                        #reset: (obsolete?)
                        w  = []
                        nr = 0

                        #radfiles.append( accfile )
                        return accfile
                            
                    else: #except:
                        pass
                    #pass
                else:  # do nothing
                    pass
            
            else: # date strange
                print 'strange, obdate, fcdate: ',hdate,gdate
                wa      = get_h5data(afile) 
                accfile = tofile(gfile,gdate,wa,acctime.seconds/3.6e3)
                return accfile

        #last bit, check...
        if acc:
            try:
                if w not in [[],None]:    

                    print 'eh ',afile,gfile,hdate,gdate
                    accfile = tofile(gfile,gdate,w/nr,int(acctime.seconds/3.6e3))
                    return accfile
            except:
                print 'Not in accumulated file: ',gfile
                
        print('strange, obdate, fcdate: ',hdate,gdate)
        return None #accfile #(Done?)

                    
    acctime = dt.timedelta(hours=accstep)

    dtg0,lt0 = find_date.datetest(gribfiles[0],pattern='fc')
    dtg1,lt1 = find_date.datetest(gribfiles[1],pattern='fc')
    deltime = (dtg1 + lt1) - (dtg0 + lt0)
    if deltime > acctime: acctime = deltime
    #sdate = find_date.datetest(h5files[0],pattern='hdf')[0] # start at first file
    #sdate = dt.datetime(1900,1,1,12,0) # arbitrary date in far history

    radfiles = set() #[]
    for gfile in sorted(gribfiles):

        alreadythere = 0 # check below if accumulated observations have been made already

        dtg,ltg = find_date.datetest(gfile,pattern='fc')
        gdate   = dtg + ltg
        #print 'In accumulate radar',os.path.split(gfile)[1], dtg,ltg,gdate

        #check if corresponding file already exists!
        resampled_files = glob.glob(os.path.join(outdir,'rad*grib'))
        for rf in sorted(resampled_files):
            dtr,ltr = find_date.datetest(rf,pattern='fc')
            if dtr == None: 
                print 'not a valid filename ',os.path.split(rf)[1]
                alreadythere = 1
                break
            elif dtr + ltr == gdate: 
                #print 'radar file available'
                alreadythere = 1
                radfiles.add(rf)
                break

        if not alreadythere:
            accfile = accumulate(gfile,h5files)
            print accfile
            radfiles.add(accfile)


    return radfiles

def resample_obs(h5files,gribfiles,outdir = './'): #,accstep=1,acc = True, start = 0):

    '''Finds SAF MSG data in a hdf5 file, 
    resamples it to fc grid and puts it in a grib message  '''

    if h5files == [] or gribfiles == []:
        print 'resample_msg: Files not found'
        return 0

    import find_date,gribtools,config

    # make dir to output files
    config.mkdir_p(outdir)

    # get necessary lats,lons ONCE:
    glats,glons = gribtools.get_latlon(gribfiles[0])
    rlats,rlons = init_MSG_h5data(h5files[0])

    def tofile(gfile,date,data,param=71,level=0):

        '''take grib file and hdf5 data, resample and insert into new grib file'''

        resampled_data = resample_radar(data,rlats,rlons,glats,glons)

        dtg     = dt.datetime.strftime(date,'%Y%m%d%H')
        newfile = os.path.join(outdir,'msg_'+str(dtg)+'+000.grib')
        gribtools.insert_grib(gfile, newfile, resampled_data,parameter = param, level = level, force = True)
        #gribtools.insert_grib(gfile, newfile, resampled_data,parameter = 0, level = level)

        return newfile

    msgfiles = set() #[]
    for gfile in gribfiles:

        alreadythere = 0 # check below if accumulated observations have been made already

        dtg,ltg = find_date.datetest(gfile,pattern='fc')
        gdate   = dtg + ltg
        print 'In resample_msg',os.path.split(gfile)[1], dtg,ltg,gdate

        #check if corresponding file already exists!
        resampled_files = glob.glob(os.path.join(outdir,'msg*grib'))
        for rf in resampled_files:
            dtr,ltr = find_date.datetest(rf,pattern='fc')
            if dtr == None: 
                print 'not a valid filename ',os.path.split(rf)[1]
                alreadythere = 1
                break
            elif dtr + ltr == gdate: 
                print 'MSG file available',rf
                alreadythere = 1
                #radfiles.append(rf)
                msgfiles.add(rf)
                break

        if not alreadythere: # do matching, resampling, write to file

            for afile in h5files:

                hdate = find_date.datetest(afile,pattern='hdf')[0]
            
                if hdate == gdate: # unlike radar, just take instantaneous cloud mask
                    wa      = get_h5data(afile,h5source = {'model':'MSG'}) 
                    msgfile = tofile(gfile,gdate,wa)

                    #print msgfile
                    msgfiles.add(msgfile)


    return msgfiles

def resample_msg(h5files,gribfiles,outdir = './'): #,accstep=1,acc = True, start = 0):

    '''Finds SAF MSG data in a hdf5 file, 
    resamples it to fc grid and puts it in a grib message  '''

    if h5files == [] or gribfiles == []:
        print 'resample_msg: Files not found'
        return 0

    import find_date,gribtools,config

    # make dir to output files
    config.mkdir_p(outdir)

    # get necessary lats,lons ONCE:
    glats,glons = gribtools.get_latlon(gribfiles[0])
    rlats,rlons = init_MSG_h5data(h5files[0])

    def tofile(gfile,date,data,param=71,level=0):

        '''take grib file and hdf5 data, resample and insert into new grib file'''

        resampled_data = resample_radar(data,rlats,rlons,glats,glons)

        dtg     = dt.datetime.strftime(date,'%Y%m%d%H')
        newfile = os.path.join(outdir,'msg_'+str(dtg)+'+000.grib')
        gribtools.insert_grib(gfile, newfile, resampled_data,parameter = param, level = level, force = True)
        #gribtools.insert_grib(gfile, newfile, resampled_data,parameter = 0, level = level)

        return newfile

    msgfiles = set() #[]
    for gfile in gribfiles:

        alreadythere = 0 # check below if accumulated observations have been made already

        dtg,ltg = find_date.datetest(gfile,pattern='fc')
        gdate   = dtg + ltg
        print 'In resample_msg',os.path.split(gfile)[1], dtg,ltg,gdate

        #check if corresponding file already exists!
        resampled_files = glob.glob(os.path.join(outdir,'msg*grib'))
        for rf in resampled_files:
            dtr,ltr = find_date.datetest(rf,pattern='fc')
            if dtr == None: 
                print 'not a valid filename ',os.path.split(rf)[1]
                alreadythere = 1
                break
            elif dtr + ltr == gdate: 
                print 'MSG file available',rf
                alreadythere = 1
                #radfiles.append(rf)
                msgfiles.add(rf)
                break

        if not alreadythere: # do matching, resampling, write to file

            for afile in h5files:

                hdate = find_date.datetest(afile,pattern='hdf')[0]
            
                if hdate == gdate: # unlike radar, just take instantaneous cloud mask
                    wa      = get_h5data(afile,h5source = {'model':'MSG'}) 

                    # sometimes file is corrupt or something
                    if wa == None: break

                    msgfile = tofile(gfile,gdate,wa)

                    print msgfile
                    msgfiles.add(msgfile)


    return msgfiles


#def radar2grib(h5files,gribfiles,outdir = './',start = 0):

    

if __name__ == '__main__':

    import pygrib
    import find_date


    h5files = glob.glob('./RAD_NL25*h5'); print 'h5: ',h5files
    gfiles  = glob.glob('./fc*grb');     print 'grib:',gfiles
    radfiles = accumulate_radar(sorted(h5files),sorted(gfiles),outdir = './test',acc = False,accstep = 3)
    sys.exit(1)

    h5files = glob.glob('./RAD_NL25*2012062121*h5'); print 'h5: ',h5files
    gfiles  = glob.glob('./acc20120621*grib');     print 'grib:',gfiles
    radfiles = accumulate_radar(sorted(h5files),sorted(gfiles),outdir = './test',acc = False,accstep = 3)

    h5files = glob.glob('./RAD_NL23*h5');          print 'h5: ',h5files
    gfiles  = glob.glob('./acc20120307*grib');     print 'grib:',gfiles
    radfiles = accumulate_radar(sorted(h5files),sorted(gfiles),outdir = './test',acc = False,accstep = 3)
    print radfiles; sys.exit(1)



    h5file = './RAD_NL23_PCP_NA_201203071400.h5'
    h5data = get_h5data(h5file)
    
    print 'retrieve: in ',find_date.datetest(h5file,pattern = 'hdf'),'\n max,min:',h5data.max(),h5data.min()

    rlats,rlons = init_h5data(h5file)
    grib_file   = 'fc2012030700+014.grb'
    grbs = pygrib.open(grib_file)
    hlats,hlons = grbs[1].latlons()
    grbs.close()

    newdata = resample_radar(h5data,rlats,rlons,hlats,hlons)
    print 'resample: in ',find_date.datetest(grib_file,pattern = 'fc'),'\n max,min:',newdata.max(),newdata.min()

    h5files = glob.glob('./RAD*h5')
    gfiles  = glob.glob('./fc*grb')
    radfiles = accumulate_radar(sorted(h5files),sorted(gfiles),outdir = './test',accstep = 2)
    print radfiles


    print 'testing "low-level" get_data from MSG satellite hdf5 files'
    h5files = glob.glob('./valkde/SAF*h5')
    for h in h5files[:3]:
        h5data = get_h5data(h,{'model':'MSG'})
        print 'hist:',numpy.histogram(h5data,bins=[0,1,2,3,4,5])

    print 'testing "high-level" resampling MSG satellite hdf5 files to harmonie forecast files'
    gribfiles = glob.glob('./valkde/ex*grib')
    msgfiles = resample_msg(h5files,gribfiles,outdir = './valkde')
    for mf in msgfiles: print mf
    print 'Done'
