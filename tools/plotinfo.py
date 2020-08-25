#! /usr/bin/env python

import os,sys

import numpy as np
import datetime as dt
import pygrib
#import h5py

    
def check_type(obj,key='values'):
    
    '''Check if object is a plotinfo instance, eg with attribute 'values' '''
    
    if hasattr(obj,'__dict__'):
        if obj.__dict__.has_key(key):
            return True
        else:
            return False
    else:
        return False

def read(datafile,**kwargs):
    
    my_inst = plotinfo()

    # test filetype, assume it is grib if there is no extension:
    if os.path.splitext(datafile)[1] in ('.h5','.hdf5'):
        my_inst.read_from_h5(datafile,**kwargs)
    elif os.path.splitext(datafile)[1] in ('.grb','.grib',''):
        my_inst.read_from_grib(datafile,**kwargs)

    return my_inst

class plotinfo():

    '''
    A class to hold essential (?) or possible information/metadata relating to a grib field
    to be able to plot it without having to keep track of all metadata yourself

    usage:
    initialise: example = plotinfo()
    read data:  example.read_from_grib(file,parameter,level, etc..,case='distinct_name')
    resample data on different grid, add, subtract other fields (on the same grid)
    and:.. plot!
    '''

    def __init__(self,f = None,par=61,level=0,leveltype='sfc',TR=0,model='Harmonie'):

        self.model = model

        self.grb       = None
        self.parameter = None
        self.leveltype = None
        self.level     = None
        self.TR        = None
        self.name      = None
        self.case      = None
        self.values    = None
        self.latlons   = None
        self.date      = None
        self.time      = None
        self.leadtime  = None

        self.colormap  = None
        self.colors    = None
        self.plot_range = None
        self.draw_line = None
        self.line_range = None
        self.plot_levels = None
        
        if f is not None:

            # possible grib file extensions: 
            gcol = ('grb','grib','GB')

            if os.path.splitext(f)[1] == '.h5':
                self.read_from_h5(self,f)
            elif [g in f[-8:] for g in gcol].count(True) > 0:
                self.read_from_grib(f,par=61,level=0,leveltype='sfc',TR=0)
            else:
                print 'Cannot automatically establish filetype ',f,'\ntry read_from_grib or read_from_hdf5'

    def copy(self):
        # make a (deep) copy
        cp = plotinfo()
        for attr in self.__dict__.keys():
            cp.__dict__[attr] = self.__dict__[attr]
        return cp

    def to_plotinfo_obj(self,field):
        # check if other is a plotinfo instance or just an array:
        if check_type(field):
            return field
        elif type(field) == np.ndarray and field.shape == self.values.shape:
            other = plotinfo()
            other.values = field
            return other
        else:
            print 'Field not compatible with current field'
            return None

    ## some methods to make it possible to compute using class instances:
    ## eg add solid to liquid precipitation, or compute wind from u,v etc

    def __add__(self,other):

        # check if other is really something:
        if other == None:
            return self.copy()
        elif other.__dict__.has_key('values') and other.values == []:
            return self.copy()

        other = self.to_plotinfo_obj(other)
        
        # check if there is already a value field (init by adding)
        if self.values == None: 
            self.values = other.values
            self = other.copy() #?
        # check if dimensions are the same:
        elif self.values.shape == other.values.shape:
            out = self.copy()
            out.values = self.values + other.values
            return out
        else:
            print 'Try to resample first:',self.values.shape,other.values.shape,self.values.shape == other.values.shape
            return None

    def __sub__(self,other):
        
        # check if other is really something:
        if self.values == []:
            return self.copy()
        elif other == None:
            return self.copy()
        elif type(other) in (int, float, complex):
            out = self.copy()
            out.values = self.values - other
            return out
        elif other.__dict__.has_key('values') and other.values == []:
            return self.copy()


        other = self.to_plotinfo_obj(other)

        # check if dimensions are the same:
        if self.values.shape == other.values.shape:
            out = self.copy()
            out.values = self.values - other.values
            return out
        else:
            print 'Try to resample first:',self.values.shape,other.values.shape,self.values.shape == other.values.shape
            return None

    def __mul__(self,other):
        
        # check if self is really something:
        if self.values == []:
            return self.copy()

        if type(other) in (float,int):
            out = self.copy()
            out.values = self.values * other
            return out

        # check if dimensions are the same:
        elif self.values.shape == other.values.shape:
            out = self.copy()
            out.values = self.values * other.values
            return out
        else:
            print 'Try to resample first:',self.values.shape,other.values.shape,self.values.shape == other.values.shape
            return None

    def __rmul__(self,other):
        
        # check if self is really something:
        if self.values == []:
            return self.copy()

        if type(other) in (float,int):
            out = self.copy()
            out.values = self.values * other
            return out

        # check if dimensions are the same:
        elif self.values.shape == other.values.shape:
            out = self.copy()
            out.values = self.values * other.values
            return out
        else:
            print 'Try to resample first:',self.values.shape,other.values.shape,self.values.shape == other.values.shape
            return None

    def __pow__(self,other):
        
        # check if self is really something:
        if self.values == []:
            return self.copy()

        if type(other) in (float,int):
            out = self.copy()
            out.values = np.power(self.values,other)
            return out

        # check if dimensions are the same:
        elif self.values.shape == other.values.shape:
            out = self.copy()
            print 'You sure? array to power of other array?'; return out

            out.values = np.exp(self.values,other.values)
            return out
        else:
            print 'Try to resample first:',self.values.shape,other.values.shape,self.values.shape == other.values.shape
            return None


    def resample(self,destgrid,nodata = 65535):


        ## if destgrid is an object with key "latlons", assume it is a plotinfo object:
        if hasattr(destgrid,'__dict__'):
            if destgrid.__dict__.has_key('latlons'):
                newlats,newlons = destgrid.latlons
        else: # assume it is a tuple with lats,lons
            newlats,newlons = destgrid

        import pyresample

        try:
            oldlats,oldlons = self.latlons
        except TypeError:
            print 'No plotting info to resample'
            return None

        old_def  = pyresample.geometry.GridDefinition(lons=oldlons,lats=oldlats)
        new_def  = pyresample.geometry.GridDefinition(lons=newlons,lats=newlats)

        # define image container
        #inradius   = 5000
        inradius   = 20000
        msg_con_nn = pyresample.image.ImageContainerNearest(self.values, old_def, 
                                                            radius_of_influence=inradius, fill_value = nodata)
        # project to new lats,lons
        area_con_nn    = msg_con_nn.resample(new_def)

        # assign to self:
        try:
            self.model = destgrid.model
        except:
            print 'only lats,lons, not model'

        self.values  = area_con_nn.image_data
        self.latlons = (newlats,newlons)
        
    def rr_from_reflec(self,nodata = 65535):

        # calculate rainrate from reflectivity
        
        data = self.values #; print data
        nodata_crit = [data[::-1,:]<255,data[::-1,:] == 255]
        nodata_do   = [np.power(10,(data[::-1,:] - 109.)/32.),nodata]

        self.values =  np.select(nodata_crit,nodata_do)
        #print self.values

    def rr_from_rac(self,nodata = 65535):

        # distill rainrate from 3h accumulated files
        
        data = self.values 
        print  'plotinfo: rr_from_rac, before: ',data.min(), data[data < nodata].mean(), data[data < nodata].max()

        nodata_crit = [data[::-1,:]<nodata,data[::-1,:] == nodata]
        nodata_do   = [0.01 * data[::-1,:],nodata]

        self.values =  np.select(nodata_crit,nodata_do)
        #print  'plotinfo: rr_from_rac, after: ',self.values.min(), self.values[self.values < nodata].mean(), self.values[self.values < nodata].max()

    def cloudmask_from_MSG(self,nodata=65535):

        data = self.values

        wn = np.select([data<1,data == 1,data==2,data==3,data>3],[nodata,0,0.1,1.,nodata])
        self.values = wn[::1,:] 

        return 0

    def CTT_from_MSG(self,nodata=65535):

        data = self.values

        '''
        calculate cloud top temperature from radiance: from computation in model, K * 0.01 in MSG!
        '''
        #wn = np.select([data<1,data>=1],[nodata,np.power(-1 * 17636684 * (data/100.), 1/4.)  - 273.15])
        wn = np.select([data<1,data>=1],[nodata,data/100.  - 273.15])
        self.values = wn[::1,:] 

        return data


    def generate_radar_latlons(self): #(dim0,dim1,ll_lon,ll_lat,ur_lon,ur_lat):

        '''
        Generate lats,lons for polar stereographically 
        projected radar from data characteristics 
        (predetermined cases: new Dutch radar, old Dutch radar, European composite)
        '''

        import pickle

        dim0,dim1 = self.values.shape

        if   dim0 == 765: ## high-res dutch radar
            raddomain = 'nlrad'
            (ll_lat,ll_lon) = (49.362,  0. )
            (ur_lat,ur_lon) = (55.389, 10.856)
        elif dim0 == 256: ## low-res dutch radar
            raddomain = 'oldrad'
            (ll_lat,ll_lon) = (49.769,  0.)
            (ur_lat,ur_lon) = (54.818,  9.743)
        elif dim0 == 512: ## european composite
            raddomain = 'eurrad'
            (ll_lat,ll_lon) = (41.937, -9.271)
            (ur_lat,ur_lon) = (58.089, 20.453)

        home = os.getenv('HOME')    
        mapdir = os.path.join(home,'python/tools')
        mappkl = os.path.join(mapdir,raddomain+'.pkl')
        if os.path.exists(mappkl):
            f = open(mappkl,'r')
            polmap = pickle.load(f)
            f.close()
            
        elif 0: # possibly faster (!?)
            sys.path.append(mapdir) #'/usr/people/plas/python/tools'
            import bmap

            polmap    = bmap.myMap(domain = raddomain,modelname = 'radar').bmap
 
        else: # if this fails (?)
            # import basemap:
            from mpl_toolkits.basemap import Basemap # map functionality

            polmap = Basemap(projection = 'stere',
                           lat_0 = 90.,  lon_0 = 0.,
                           lat_ts = 60., 
                           llcrnrlon = ll_lon, llcrnrlat = ll_lat,
                           urcrnrlon = ur_lon, urcrnrlat = ur_lat 
                           )
        
        lons,lats = polmap.makegrid(dim1,dim0,returnxy=False)
        #print 'LON',lons,'\nLAT',lats
        self.latlons = (lats,lons)

    def generate_sat_latlons(self,saf=True,verb=False):
        
        '''
        Generate lats,lons for MSG geostationary satellite image from map characteristics 
        Possible to read some stuff from the h5 file, not necessary for now...
        '''

        dim0,dim1 = self.values.shape
        print dim0,dim1

        import pyproj

        projs = '+proj=geos +a=6378.1690 +b=6356.5838 +lat_0=0.0 =+lon_0=0.0 +h=35785.8310'
        p = pyproj.Proj(projs)

        if (dim1,dim0) == (3712,928): #saf:
            x_lr, y_lr =  5565.7482275008615, 2655.3569710718393
            x_ul, y_ul = -5568.748630858005,  5436.730883143699
            
        elif (dim1,dim0) == (3000,1000): #METEOSAT
            x_lr, y_lr =  4500.5805, 2400.3096 
            x_ul, y_ul = -4500.5805, 5400.6966 
            pix = 3.00387 #?

        xs = np.linspace(x_ul, x_lr, dim1)
        ys = np.linspace(y_ul, y_lr, dim0)
        xn,yn = np.meshgrid(xs,ys)

        lons,lats = p(xn,yn,inverse=True)
        print 'Generated lons, lats: ',lons.shape, lats.shape, self.values.shape
        self.latlons = (lats,lons)

        return 0
 

    def generate_sat_latlons_old(self,saf=True,verb=False):

        '''
        Generate lats,lons for MSG geostationary satellite image from map characteristics 
        recognize by keyword SAF, which has a slightly different map
        '''
        
        # import basemap:
        from mpl_toolkits.basemap import Basemap # map functionality

        dim0,dim1 = self.values.shape

        # from hdf5 file:
        # proj string: +proj=geos +a=6378169.0 +b=6356583.8 +lon_0=0.0 +h=35785831.0

        if saf: # from hdf5 file
            x_lr =  5565748.2275008615
            x_ul = -5568748.630858005
            y_lr =  2655356.9710718393
            y_ul =  5436730.883143699

        else:
            ## msg METEOSAT_10 proj:
            '''
            geo_column_offset = -1500.0
            geo_row_offset = -1800.0
            geo_pixel_size_x,geo_pixel_size_y = 3.000387,3.000387
            
            x = (3000 + geo_column_offset) * geo_pixel_size_x
            y = (1000 + geo_row_offset) * geo_pixel_size_y
            '''
        
            x_lr =   4500.5805
            x_ul =  -4500.5805 
            y_lr =   2400.3096
            y_ul =   5400.6966
            projs = '+proj=geos +a=6378.1690 +b=6356.5838 +lat_0=0.0 =+lon_0=0.0 +h=35785.8310'

            pix = 3.00387

            xs = np.linspace(x_ul + pix/2, x_lr-pix/2, 3000)
            ys = np.linspace(y_ul+pix/2, y_lr-pix/2, 1000)
            xn,yn = np.meshgrid(xs,ys)

            import pyproj
            p = pyproj.Proj(projs)
            lons,lats = p(xn,yn,inverse=True)
            self.latlons = (lats,lons)
            return 0

        maph = Basemap(projection = 'geos',
                       rsphere=(6378169.00,6356583.8),
                       lon_0=0.0,
                       llcrnrx = x_ul,llcrnry = y_lr,
                       urcrnrx = x_lr,urcrnry = y_ul,
                       satellite_height=35785831.0
                       )

        xtot = x_lr - x_ul
        ytot = y_ul - y_lr
        xs = np.arange(x_ul + xtot/(2*dim1),x_lr,xtot/dim1)
        ys = np.arange(y_lr + ytot/(2*dim0),y_ul,ytot/dim0)

        if verb: print len(xs),len(ys),min(xs),x_ul,min(ys),max(ys) #; sys.exit(1)
        
        x  = np.array([xs for i in range(dim0)])
        y  = np.array([ys for i in range(dim1)]).transpose()

        print 'Image shape ', x.shape, y.shape

        lonpt, latpt  = maph(x,y,inverse=True); #print lonpt,latpt
        #lons,lats,x,y = maph.makegrid(dim1,dim1,returnxy=True) #False
        lons,lats,x,y = maph.makegrid(dim0,dim1,returnxy=True) #False
        
        if verb:
            print '4kant: ',lons.shape, lons[100,1856],lons[-100,1856], lons[1856,100],lons[1856,-100]
            print '4kant: ',lats.shape, lats[100,1856],lats[-100,1856], lats[1856,100],lats[1856,-100]#; sys.exit(1)
            print lonpt.shape, latpt.shape

            a,b = 1856,1856
            print 'zentrum       ', lons[a,b], lats[a,b],x[a,b],y[a,b]

        #return lats[-928:,:],lons[-928:,:]
        #self.latlons = (lats[-928:,:],lons[-928:,:])
        self.latlons = (lats[:,:],lons[:,:])
        return 0

    def generate_rawsat_latlons(self,saf=True,verb=True):

        '''Generate lats,lons for MSG geostationary satellite image from map characteristics '''
        
        # import basemap:
        from mpl_toolkits.basemap import Basemap # map functionality

        dim0,dim1 = self.values.shape

        # from hdf5 file:
        # proj string: +proj=geos +a=6378169.0 +b=6356583.8 +lon_0=0.0 +h=35785831.0

        if saf: # from hdf5 file
            x_lr =  5565748.2275008615
            x_ul = -5568748.630858005
            y_lr =  2655356.9710718393
            y_ul =  5436730.883143699

        else:
            ## msg METEOSAT_10 proj:
            '''
            geo_column_offset = -1500.0
            geo_row_offset = -1800.0
            geo_pixel_size_x,geo_pixel_size_y = 3.000387,3.000387
            
            x = (3000 + geo_column_offset) * geo_pixel_size_x
            y = (1000 + geo_row_offset) * geo_pixel_size_y
            '''
        
            x_lr =   4500580.5
            x_ul =  -4500580.5 
            y_lr =   2400309.6
            y_ul =   5400696.6
            


        maph = Basemap(projection = 'geos',
                       rsphere=(6378169.00,6356583.8),
                       lon_0=0.0,
                       llcrnrx = x_ul,llcrnry = y_lr,
                       urcrnrx = x_lr,urcrnry = y_ul,
                       satellite_height=35785831.0
                       )

        xs = np.linspace(x_ul, x_lr, 3001)
        ys = np.linspace(y_ul, y_lr, 1001)
        xn,yn = np.meshgrid(xs,ys)
        lonpt, latpt  = maph(xn,yn,inverse=True);

        xtot = x_lr - x_ul
        ytot = y_ul - y_lr
        xs = np.arange(x_ul + xtot/(2*dim1),x_lr,xtot/dim1)
        ys = np.arange(y_lr + ytot/(2*dim0),y_ul,ytot/dim0)

        if verb: print len(xs),len(ys),min(xs),x_ul,min(ys),max(ys) #; sys.exit(1)
        
        x  = np.array([xs for i in range(dim0)])
        y  = np.array([ys for i in range(dim1)]).transpose()

        print 'Image shape ', x.shape, y.shape

        lonpt, latpt  = maph(x,y,inverse=True); #print lonpt,latpt
        #lons,lats,x,y = maph.makegrid(dim1,dim1,returnxy=True) #False
        lons,lats,x,y = maph.makegrid(dim1,dim0,returnxy=True) #False

        #self.latlons = (lats[-928:,:],lons[-928:,:])
        self.latlons = (lats[:,:],lons[:,:])
        return self.latlons


    def get_settings(self,preset = None):

        '''
        an interface to predefined settings: colors,levels,Celsius/knots etc
        '''

        if os.path.isdir('/nobackup_1/users/plas'):
            fabriek = '/nobackup_1/users/plas/python/fabriek/fabriek'
            fabriek = '/usr/people/plas/python/fabriek'
        elif os.path.isdir('/home/plas/fabriek'):
            fabriek = '/home/plas/fabriek'
        else:
            print 'find settings file'
            fabriek = '/path/to/settingsfile'
        sys.path.append(fabriek)
        import settings

        # allow for sloppy settings:
        if preset.lower() in [a.lower() for a in ('PCP','APCP','NCPCP','precip','RR','pcpacc','pcp_acc','pcpa')]:
            preset = settings.PCP_acc
        elif preset.lower() in [a.lower() for a in ('PCPI','precip_intensity','pcp_i','precip_int','RRI')]:
            preset = settings.PCP_i
        elif preset.lower() in ('t2m','t','tsurf','temp2m','temp','temperature'):
            preset = settings.T2M
            if self.values.max() > 100: self.values -= 273.15 
        elif preset in ('U10','V10','uv','UV','wind'):
            preset = settings.WIND
        elif preset in ('P','MSLP','pres','PRES','press','PRESS','pressure'):
            preset = settings.PRESS
        elif preset in ('CLC','clc','tcc','cloud','CLOUD','cloudcover'):
            preset = settings.CLOUD
        elif preset in ('VIS','vis','mist','soep','zicht'):
            preset = settings.VIS
        elif preset == None:
            preset = settings.default

        #if not self.colormap = None
        print preset
        for k,v in preset.iteritems():
            print k,v
            self.__dict__[k] = v
        

    def grib2class(self,grb,name=None):

        '''
        strips some metadata from the header, adds it to instance
        '''

        self.grb       = grb
        self.parameter = grb.indicatorOfParameter #par
        self.leveltype = grb.indicatorOfTypeOfLevel #leveltype
        self.level     = grb.level
        self.TR        = grb.timeRangeIndicator #TR
        if name is not None: 
            self.name  = name
        else: 
            self.name = grb.name
        self.values   = grb.values; #print 'in plotinfo: empty?, param, level,case:',self.values == [], par,level,case #self.values
        self.latlons  = grb.latlons()
        self.date     = grb.dataDate
        self.time     = grb.dataTime/100
        self.leadtime = grb.endStep # in hours...
        
        self.datetime = dt.datetime.strptime(str(self.date), '%Y%m%d') + dt.timedelta(hours=self.time) + dt.timedelta(hours=self.leadtime)


    def read_from_grib(self,gribfile,par = 61,level = 0,leveltype='sfc',TR=0,name=None,case='test',model = 'Harmonie'):

        '''
        Read a grib file and derive parameters necessary to plot (or inspect) the file
        '''
        

        ## instantiate
        #self.__init__(model)
        self.values = []

        self.model    = model
        self.case     = case

        # possibly translate leveltype
        ltrans = {#105:'heightAboveGround',
                  105:'sfc',
                  #109:'hybrid',
                  109:'ml',
                  102:'heightAboveSeaLevel'
                  }

        if type(leveltype) == int and leveltype in ltrans.keys():
            leveltype = ltrans[leveltype]

        #import pygrib
        if type(gribfile) == pygrib.index:
            try: # based on all four criteria:
                # if 1:
                selected_grbs=gribfile.select(
                    indicatorOfParameter   = par,
                    indicatorOfTypeOfLevel = leveltype,
                    level                  = level,
                    timeRangeIndicator     = TR
                    )

            except: # if leveltype is "different"

                print 'Not found: ',par,leveltype,level,TR 
                return None       
 
            grb = selected_grbs[0]
            self.grib2class(grb)

        else:
            grbs = pygrib.open(gribfile)
            for grb in grbs:
                if grb.indicatorOfParameter == par \
                        and grb.levelType == leveltype \
                        and grb.level == level \
                        and grb.timeRangeIndicator == TR:
                
                    print 'found!',grb
                    self.grib2class(grb)

                    break # If one is found, then break out of loop

            # clean up
            grbs.close()

        #print 'Time ?', self.date, self.time
        if self.values == []: print 'Not found',par,level

        return self # experiment


    def read_from_h5(self,h5file,
                     datapath = 'image1/image_data',metapath='overview',
                     dtattr='image_datetime_valid',
                     name='precipitation',model='obs'):

        
        import h5py
        import datetime as dt

        print 'checking path' , name.lower(),'msg' in name.lower(),'ir' in name.lower()
        if 'msg' in name.lower():
            isMSG   = True
            isMSGCM = False
            isMSGIR = False
            if 'ir' in name.lower():
                isMSGIR = True
                print isMSGIR
            elif 'cm' in name.lower():
                isMSGCM = True
            else:
                print 'unknown MSG product',name.lower()
        else:
            isMSG   = False
            isMSGIR = False
            isMSGCM = False

        try:
            f=h5py.File(h5file,'r')
        except IOError:
            print 'File corrupt: ',h5file
            return 0

        try:
            # once you know where the data is, fill in 'path':
            data = f[datapath]
            self.values = data[:,:] #/100
            (dim0,dim1) = self.values.shape; print 'found data ',dim0,dim1
        except IOError:
            print 'Reading file failed',h5file,datapath
            return 0

        # get valid date,time from metadata path
        # and convert to datetime object and then to int
        if isMSG:
            if isMSGCM:
                dts = f.attrs[u'TIME_STAMP_UP_LINE'] # eg '20131015180908'
                dtd = dt.datetime.strptime(dts,'%Y%m%d%H%M%S') # eg '21-JUN-2012;12:00:00.000'
                print dtd
            elif isMSGIR: # eg image_acquisition_time = 31-JAN-2014;03:00:00.000
                dts = f[metapath].attrs[u'image_acquisition_time'][0] 
                print dts, type(dts)
                dtd = dt.datetime.strptime(dts,'%d-%b-%Y;%H:%M:%S.000') 
                print dtd
        else: # Radar
            print 'radar:',metapath,dtattr
            try:
                dts = f[metapath].attrs[dtattr][0] #['product_datetime_start'][0]
                dtd = dt.datetime.strptime(dts,'%d-%b-%Y;%H:%M:%S.000') # eg '21-JUN-2012;12:00:00.000'
            except:
                print 'Problem:',f[metapath]
                dtd = dt.datetime(2004,5,14,16)

        date = dtd.strftime('%Y%m%d')
        time = dtd.strftime('%H')

        # clean up
        f.close()

        self.case = model

        # some custom methods, attributes, use decorator?
        if isMSGCM:
            self.cloudmask_from_MSG()
            self.generate_sat_latlons(saf=True)
            #self.case = 'MSG'
            self.parameter = 71
            self.colormap = 'binary_r'
        elif isMSGIR:
            dd = self.CTT_from_MSG()
            ll = self.generate_sat_latlons(saf=False)
            #self.case = 'MSG'
            self.parameter = 114
            self.levels   = np.linspace(-10,-40,31)
            self.colormap = 'binary'

        else:
            if 'RAC' in os.path.split(h5file)[1]:
                self.rr_from_rac()
            else:
                # compute rain rate from reflectivity, get lats,lons
                self.rr_from_reflec()
                
            # get the lats,lons from the polar stereographic projection
            self.generate_radar_latlons()
            #self.case = 'radar'
            self.parameter = 61


        self.level = 0
        self.leveltype = 'sfc'

        #self.values = self.rr_from_reflec(w)
        self.datetime = dtd
        self.date   = date
        self.time   = time
        self.leadtime = 0
        self.name   = name
        self.nodata = 65535
        
        return 0 #dd,ll  #pass

    def plot(self,out=None,domain = 'nl',lsmask_colour = 'black'):


        if self.values  == []: return None
        if self.latlons == None: return None

        sampledir = '/usr/people/plas/python/tools'

        if 1: # to do: test if domain is not some other custom definition
            if domain == 'nl':
                grbs = pygrib.open(os.path.join(sampledir,'samplegrib_SA.grb'))
                nlgrid = grbs[1].latlons()
                #print nlgrid; sys.exit(0)
                self.resample(nlgrid)
                grbs.close()
            elif domain == 'eur':
                grbs = pygrib.open(os.path.join(sampledir,'samplegrib_LA.grb'))
                eurgrid = grbs[1].latlons()
                self.resample(eurgrid)
                grbs.close()
        else: #
            print "resample failed, try original grid: "



        import matplotlib.pyplot as plt
        from matplotlib import cm as CM
        from matplotlib import rcParams
        rcParams['contour.negative_linestyle'] = 'solid' # solid lines for negative contours

        home = os.getenv('HOME')
        sys.path.append(os.path.join(home,'python/tools'))
        import bmap
        #print bmap.__file__; sys.exit(0)
        
        # create figure canvas, axes
        fig  = plt.figure()
        ax   = fig.add_subplot(1,1,1)

        # create mapping, draw land-sea mask, meridians etc
        mplot = bmap.myMap(domain=domain,dbg = True)
        mplot.xymap(latlons = self.latlons, domain = domain) #, modelname = self.model)
        #print 'x,y map:',domain,mplot.x.shape, mplot.y.shape
        mplot.set_axes(ax)
        mplot.dress_up_map(domain = domain,lsmask_colour = lsmask_colour)

        #mplot.bmap.fillcontinents(color='Wheat',lake_color='#99ffff')
        #mplot.bmap.drawmapboundary(fill_color='aqua')


        
        exarg = 'max' #'both' # extending colors beyond the given levels, could be 'both, 'neither', see matplotlib doc

        #### where to put this kind of intelligence?
        if self.plot_levels == None:
            if self.plot_range: # is not None: make range of plot levels
                r = self.plot_range # for brevity:
                try:
                    self.plot_levels = np.arange(r[0],r[1],r[2])
                except: 
                    print 'Not a valid range definition: ',r #,'\nUsing default (0,10,1)'
            else:
                gmax,gmin = self.values.max(),self.values.min()
                self.plot_levels = np.arange(gmin,gmax,(gmax-gmin)/10.)
                print gmin,gmax,self.plot_levels

        if self.colormap == None and self.colors == None: 
            self.colormap = CM.get_cmap('jet')
        elif self.colors is not None:
            self.colormap = None
        elif self.colormap is not None:
            cmap = self.colormap
            try:
                self.colormap = CM.get_cmap(cmap)
            except: 
                print 'Is not available? See matplotlib doc ',self.colormap
                sys.exit(0)
            #cmap.set_over(color='k',alpha = None) # what, why?
            self.colors = None
        else:
            pass

        if self.__dict__.has_key('lines'): 
            lines = self.lines
            print lines
        elif self.__dict__.has_key('line_range'):
            if self.line_range == None: pass
            r = self.line_range; print r
            try:
                lines = np.arange(r[0],r[1],r[2])
            except:
                print 'Not a valid range definition: ',r
        else:
            print 'contourplot draw_line: using levels ',self.plot_levels
            lines = self.plot_levels # probably redundant
            
        # check if name is well-behaved:
        # characters not allowed in strings for file names:
        import re
        re_ill =  re.compile(r"^[^<>/{}[\]~`\s\*\\\?\|]*$");

        if re_ill.match(str(self.case)):
            #print("RE2: All chars are valid.",self.case)
            pass
        else:
            print("RE2: Not all chars are valid.",self.case)
            sys.exit(1)
 
        ##########

        if 1:
            print 'shapes: values, latlons, xy', self.values.shape, self.latlons[0].shape, self.latlons[1].shape, mplot.x.shape, mplot.y.shape

        pcont = mplot.bmap.contourf(mplot.x,mplot.y,self.values,
                                       self.plot_levels,cmap=self.colormap,colors = self.colors,extend=exarg,ax=ax)# ,zorder=100) # for filled continents

        if self.draw_line: 
           plcont = mplot.bmap.contour(mplot.x,mplot.y,self.values,lines,colors='black',ax=ax) #,zorder=101) # for filled continents


        # create colorbar
        cb = fig.colorbar(pcont,shrink=0.9, extend=exarg,format='%.3g') #format='%.1e') 
        cb.set_label(self.name,fontsize=9)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(9)

        # title
        print 'setting title',self.date
        title  = '{mod}: par {parltr} \n at {date}+{lead}, validtime {validt}'
        ax.set_title(title.format(mod= self.model, 
                                  parltr = str(self.parameter)+':'+str(self.leveltype)+':'+str(self.level), 
                                  date = str(self.date)+str(self.time).zfill(2),lead = str(self.leadtime).zfill(3),
                                  validt = self.datetime.strftime('%Y%m%d%H')) )


        # save plot
        print 'setting filename',self.date
        #stdoutf = './plot_'+str(self.case)+'_'+str(self.parameter)+'_'+str(self.leveltype)+'_'+str(self.level)+'_'+str(self.date)+str(self.time).zfill(2)+'+'+str(self.leadtime).zfill(3)+'_'+domain+'.png'
        stdoutf = './plot_{case}_{param}_{ltype}_{l}_{dtg}+{lt}_{dom}.png'.format(case=self.case,param=self.parameter,ltype=self.leveltype,l=self.level,dtg=str(self.date)+str(self.time).zfill(2),lt=str(self.leadtime).zfill(3),dom=domain)

        if out == None:
            out = stoutdf
        elif os.path.isdir(out):
            out = os.path.join(out,stdoutf)
        else:
            print 'create ',out

        fig.savefig(out, dpi = 100)

        # for many plots a memory leak...
        #del(ax)
        #fig.clf()
        #del(fig)
        plt.close(fig) ## <- this is it! This closes the plot instance.

        print 10*'=','> created ',out
        return out

    def create_histogram(self,bins=[0,0.1,0.3,1,3,5,10,15,25,50,100,150,250],outdir = './',verb=False):

        ## import stuff
        import matplotlib.pyplot as plt

        myhist,bins = np.histogram(self.values.ravel(),bins = bins)
    
        # create figure canvas, axes
        fig  = plt.figure()
        ax   = fig.add_subplot(1,1,1)  

        if verb:
            print range(len(bins)-1)
            for b,n in zip(bins,myhist): print b,n
        ax.bar(range(len(bins)-1),myhist)

        xTickMarks = [str(b) for b in bins]
        if verb: print xTickMarks
        xtickpos   = ax.set_xticks([i+0.3 for i in range(len(bins)-1)])
        xtickNames = ax.set_xticklabels(xTickMarks)

        if not hasattr(self,'ltbeg'): self.ltbeg = 0
        if not hasattr(self,'ltend'): self.ltend = self.leadtime
    
        histout = 'hist_{model}_{dtg}_{st}_{ltbeg}_{ltend}.png'.format(model=self.model,
                                                                       dtg=str(self.date)[0:6],
                                                                       st=str(self.time).zfill(2),
                                                                       ltbeg=str(self.ltbeg).zfill(3),
                                                                       ltend=str(self.ltend).zfill(3)
                                                                       )
        out = os.path.join(outdir,histout)
        fig.savefig(out, dpi = 200)
        print 10*'=','> created ',out


    def savegrib(self,newfile,newparam=None,fromgrib=None):

        
        if not hasattr(self,'grb') or self.grb == None:
            if fromgrib==None: 
                print 'Need example grib to encode to grib (for now)'; sys.exit(1)
            else:
                if os.path.exists(fromgrib):
                    mgrb = plotinfo(model='model')
                    try:
                        mgrb.read_from_grib(fromgrib,61,457)
                        print 'Found values ',mgrb.values.shape
                    except:
                        mgrb.read_from_grib(fromgrib,181,0)
                        print 'Found values ',mgrb.values.shape
                    
                    if mgrb.latlons and self.latlons: # test if all lats,lons are available
                    
                        self.resample(mgrb)
                        self.grb = mgrb.grb
                        print 'Resampled values ',mgrb.values.shape
                    else:
                        print 'Problem reading lats,lons'
                else:
                    print 'Path does not exist:',fromgrib
                    sys.exit(1)

        # NB if eg nodata is too large, values in grib message will be "quantized"!!
        # make sure that really large values are not really really large
        self.values[self.values >= 1000] = 1000

        # make a new grib file
        newgrib = open(newfile,'wb')
        grbref        = self.grb
        grbref.values = self.values

        msg = grbref.tostring()
        newgrib.write(msg)

        # close files
        newgrib.close()
        print 'Wrote ',newfile


if __name__ == '__main__':

    home = os.getenv('HOME')
    datadir = '/Users/plas/HARP/Harp/data'
    datadir = os.path.join(home,'HARP/Harp_sample/data')

    # define some example 
    gribfile_prev = os.path.join(datadir,'fc','harm','ex2012062100+011grib')
    gribfile      = os.path.join(datadir,'fc','harm','ex2012062100+012grib')
    radfile       = os.path.join(datadir,'radar','eur_15m','RAD_NL23_PCP_NA_201206211200.h5')
    radfilenl     = os.path.join(datadir,'radar','nl_3h','RAD_NL25_RAC_03H_201206212000.h5')

    if 0: # debug: print out keys,values
        grbs = pygrib.open(gribfile)
        for k in grbs[1].keys():
            print k, grbs[1][k]
        sys.exit(0)

    gr11_new = read(gribfile_prev,par=61,level=457,model = 'test')
    gr11_new.get_settings(preset = 'APCP')
    gr11_new.plot(out = home,domain = 'eur',lsmask_colour='black')
    sys.exit(0)

    gr11 = plotinfo(model='test1')
    gr11.read_from_grib(gribfile_prev,61,457,model = 'test')

    print 'opening ',gribfile,os.path.exists(gribfile)
    gr12 = plotinfo(model='test1')
    gr12.read_from_grib(gribfile,61,457,model = 'test')
    gr12.get_settings(preset = 'APCP')
    print dir(gr12)

    somevalues = gr11.values
    gr12 = gr12 + somevalues

    # deduct previous timestep
    hour = gr12 - gr11
    hour.get_settings(preset = 'APCP')
    
    # test outdir/outfile functionality
    hour.plot(out = home,domain = 'eur',lsmask_colour='black')
    #sys.exit(0)

    gr12.plot(out = 'test_plinf.png')
    gr12.plot(out = 'test_eurplinf.png',domain='eur')


    print 'opening ',radfile,os.path.exists(radfile)
    rad12 = plotinfo(model='testrad')
    rad12.read_from_h5(radfile)
    rad12.get_settings(preset = 'APCP')
    print dir(rad12)
    #sys.exit(1)

    rad12.plot(out = 'test_eurradplinf.png',domain = 'eurrad',lsmask_colour = 'black')
    print 'radar rainrate values: ', rad12.values[rad12.values < 500].mean(),rad12.values[rad12.values < 500].max()

    print 'opening ',radfilenl,os.path.exists(radfilenl)
    radnl = plotinfo(model='testrad')
    radnl.read_from_h5(radfilenl)
    radnl.get_settings(preset = 'APCP')
    #print dir(radnl)
    radnl.plot(out = 'test_nlradplinf.png',domain = 'nlrad',lsmask_colour = 'black')
    print 'radar rainrate values: ', radnl.values[radnl.values < 500].mean(),radnl.values[radnl.values < 500].max()
    #sys.exit(0)

    # test resampling, subtraction
    radfc = rad12.copy()
    print 'before',radfc.values.shape
    radfc.resample(gr12)
    print 'after',radfc.values.shape
    radfc.plot('test_eurradres.png',domain = 'eur')
    #radfc.values = gr12.values - radfc.values
    diff = hour - radfc

    effdiff = diff.values[abs(diff.values) < 2000]
    MAE = abs(effdiff).sum()/effdiff.size
    print 'MAE (sum(fc - ob)/#participating grid points)',MAE

    MSE = (effdiff**2).sum()/effdiff.size
    print 'MSE  (sum(fc - ob)^2)/(#participating grid points)',MSE
    RMSE = np.sqrt((effdiff**2).sum())/effdiff.size
    print 'RMSE sqrt(sum(fc - ob)^2)/(#participating grid points)',RMSE

    if 1:
        diff.colors = None
        diff.colormap = 'RdBu_r'
        scale = abs(diff.values[abs(diff.values) < 500]).max()/20.
        print scale
        scale = 5.
        diff.plot_levels =  np.arange(-scale,scale+0.1,(scale)/10.)
        diff.plot(out = 'test_eurraddiff.png',domain = 'eur')
        # a lot of lines:

        #diff.draw_line   = True
        #diff.line_range  =  np.arange(-scale,scale+0.1,2.)
        #diff.plot(out = 'test_nlraddiff.png',domain = 'nl')
        #sys.exit(0)

    # test other parameter:
    t12 = plotinfo(model='test1')
    t12.read_from_grib(gribfile,11,2,model = 'harm')
    t12.get_settings(preset = 'T2M')
    t12.plot(out = 'test_t2mplinf.png',domain = 'eur',lsmask_colour = 'purple')

    msgfile = 'SAF_MSG3_CMa_201308060400.h5'
    msg = plotinfo(model='msg')
    msg.read_from_h5(msgfile,name = 'MSG') # NB: naming with something like msg is important to recognize the file...

    # to plot, resample to 'known' bmap grid, such as a harmonie grid
    msg.resample(gr12)
    msg.plot(domain = 'nl')

    msgirfile = 'SAF_MSG3_CMa_201308060400.h5'
    msg = plotinfo(model='msg ir')
    msg.read_from_h5(msgfile,name = 'MSGIR') # NB: naming with something like msg is important to recognize the file...

    # to plot, resample to 'known' bmap grid, such as a harmonie grid
    msg.resample(gr12)
    msg.plot(domain = 'nl')
