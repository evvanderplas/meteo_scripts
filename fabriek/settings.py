#/usr/bin/env python2.6

def init_list(model = 'Harmonie', version = '3712'):

    print 'init_list in settings: OBSOLETE!'; sys.exit(1)

    gribtable = {}
    if model == 'Harmonie' and int(version) >= 3712: 
        tabfile = 'tools/2.233.253.table'

    else:
        # use old gribtab
        gribtable['p'] = 1
        gribtable['r'] = 61
        gribtable['s'] = 62
        gribtable['gr'] = 63
        gribtable['u'] = 33
        gribtable['v'] = 34
        gribtable['t'] = 11
        gribtable['tcc'] = 71
        gribtable['clwc'] = 76
        gribtable['ciwc'] = 58
        gribtable['162'] = 162
        gribtable['163'] = 163
        
        gribtable['acclevel'] = (457,0)
        gribtable['intlevel'] = (456,0)
        gribtable['gustlevel'] = (10,0)
        return gribtable

    f= open(tabfile,'r')
    for line in f:
        s = line.split()
        gribtable[s[1]] = int(s[0])

        # not elegant:
        gribtable['acclevel'] = (0,4)
        gribtable['intlevel'] = (0,0)
        gribtable['gustlevel'] = (10,2)

    f.close()
    return gribtable

#print 'In settings:',source; sys.exit(1)
#gribtable = init_list()
#print gribtable

default = {'name':'Field',
           'shortname':'fld',
           'grib_indicator':0,
           'level_indicator': 105, #103?
           'level': 0, 
           'tr_ind':0,
           #gribtable['intlevel'][1],
           #'plot_range': (90e3,110e3,1e3), 
           #'draw_line':True,
           #'line_range':(90.e3,110.e3,5.e2)
           'colormap': 'jet'
           }
PRESS = {'name':'Pressure above groundlevel',
         'shortname':'p',
         'grib_indicator':'p', #gribtable['p'],
         'level_indicator': 105, #103?
         'level': 0, 
         'plot_range': (90e3,110e3,1e3), 
         'colormap': 'jet',
         'draw_line':True,
         'line_range':(90.e3,110.e3,5.e2)
       }
T2M = {'name':'2 metre temperature',
       'shortname':'t2m',
       'grib_indicator':11, #gribtable['t'],
       'level_indicator': 105,
       'level': 2, 
       'plot_range': (-10,30,1), #(260,310,1), 
       'plot_range': (16,41,2), #(260,310,1), 
       'colormap': 'jet',
       'draw_line':True,
       'line_range':(-0,26,5)
       }
BRTEMP = {'name':'Cloud top temperature',
          'shortname':'tcl',
          'grib_indicator':118,
          'level_indicator': 8,
          'level': 39680, 
          'grib_indicator_n':114,
          'level_indicator_n': 8,
          'level_n': 0, 
          'tr_ind_n': 0, 
          'plot_range': (-40,10.1,2), #(260,310,1), 
          'colormap': 'gray_r'
          }
RELH = {'name':'Relative humidity',
        'shortname':'rh',
        'grib_indicator':52,#gribtable['relh'],
        'level_indicator': 105,
        'level': 2, 
        'plot_range': (0,1.15,.1), 
        'colormap': 'Blues'
       }
Q2M  = {'name':'Specific humidity',
        'shortname':'q2m',
        'grib_indicator':51,#gribtable['relh'],
        'level_indicator': 105,
        'level': 2, 
        'plot_range': (0,1.0001e-2,5.e-4), 
        'colormap': 'jet_r'
       }

VV  = {'name':'vertical wind',
       'shortname':'vv',
       'grib_indicator':40, #gribtable['u'],
       'level_indicator': 109,
       'level': 35, 
       'plot_range': (-6,7,1), 
       'colormap': 'RdBu_r',
       'draw_line':True,
       #'line_range':(-10,10,5)
       'lines':(-6,-3,3,6)
       }

U10 = {'name':'10 metre meridional wind',
       'shortname':'u10',
       'grib_indicator':33, #gribtable['u'],
       'level_indicator': 105,
       'level': 10, 
       'plot_range': (.01,21,1), 
       'colormap': 'Reds'
       }
V10 = {'name':'10 metre longitudinal wind',
       'shortname':'v10',
       'grib_indicator':34, #gribtable['v'],
       'level_indicator': 105,
       'level': 10, 
       'plot_range': (.01,21,1), 
       'colormap': 'Reds'
       }
UGST = {'name':'10 metre meridional gust',
        'shortname':'ugst',
        'grib_indicator':162, #gribtable['162'],
        'level_indicator': 105,
        'level': 10, 
        'tr_ind':2,#gribtable['gustlevel'][1], # either 0 or 2
        'plot_range': (.1,21,1), 
        'colormap': 'Reds'
        }
VGST = {'name':'10 metre longitudinal gust',
        'shortname':'vgst',
        'grib_indicator':163, #gribtable['163'],
        'level_indicator': 105,
        'level': 10, 
        'tr_ind':2, #gribtable['gustlevel'][1], # either 0 or 2,
        'plot_range': (.1,21,1), 
        'colormap': 'Reds'
        }
WIND = {'name':'10 metre wind', # (kts)',
        'shortname':'wind',
        'grib_indicator':33, #[33,34], #dummy, only used for plot settings
        'level_indicator': 105,
        'level': 10, 
        #'plot_range': (.1,21,1), 
        #'colormap': 'Reds'
        #'plot_levels':(0,7,11,16,21,27,34,41,48,56,64,74,86,100,200),
        'plot_levels': (0,2, 4, 6, 8,10,12,14,16,18,20,22,24, 26, 28),
        'colors': ('LightCyan', #'Aquamarine',
                   'PaleGreen','LawnGreen','green','Aqua','RoyalBlue','Blue',
                   'magenta','orange','red','yellow','lightgrey','grey','black'),
        }

GUST = {'name':'10 metre gust',
        'shortname':'gust',
        'grib_indicator':162, #[162,163], #dummy, only used for plot settings
        'level_indicator': 105,
        'level': 10, 
        'tr_ind':2,
        #'plot_range': (1,41,1), 
        #'colormap': 'OrRd',
        'draw_line':True,
        'line_range':(20,41,5),
        #'colormap': 'Reds',
        #'plot_levels': (0,2, 4, 6, 8,10,12,14,16,18,20,22,24, 26,30,35,40),
        #'colors': ('LightCyan', #'Aquamarine',
        #           'PaleGreen','LawnGreen','green','Aqua','RoyalBlue','Blue',
        #           'magenta','orange','red','yellow','ivory','lightgrey','grey','Sienna','black'),
        #'plot_levels': (-5,5,15,30,45,60,75,90,105,120,135,150,165),
        'plot_range': (16,41,2),
        #'colors':( (255,255,255),(235,235,235),(190,190,190),(150,150,150),(100,200,255),(50,120,255),
        #           (0,66,204),(34,12,127),(223,100,255),(173,50,228),(123,0,178),(255,120,120),(255,0,0)
        #       )
        'colors':[[1.0, 1.0, 1.0],
                  [0.9215686274509803, 0.9215686274509803, 0.9215686274509803],
                  [0.7450980392156863, 0.7450980392156863, 0.7450980392156863],
                  [0.5882352941176471, 0.5882352941176471, 0.5882352941176471],
                  [0.39215686274509803, 0.7843137254901961, 1.0],
                  [0.19607843137254902, 0.47058823529411764, 1.0],
                  [0.0, 0.25882352941176473, 0.8],
                  [0.13333333333333333, 0.047058823529411764, 0.4980392156862745],
                  [0.8745098039215686, 0.39215686274509803, 1.0],
                  [0.6784313725490196, 0.19607843137254902, 0.8941176470588236],
                  [0.4823529411764706, 0.0, 0.6980392156862745],
                  [1.0, 0.47058823529411764, 0.47058823529411764],
                  [1.0, 0.0, 0.0]],

        }

PCP_i = {'name':'Precipitation intensity',
        'shortname':'pcpi',
         'grib_indicator': 61,#gribtable['r'],
         'level_indicator': 105,
         'level': 456,#gribtable['intlevel'][0], 
         'tr_ind':0,#gribtable['intlevel'][1],
         'grib_indicator_n': 181,
         'level_n': 0,
         'tr_ind_n':0,
         #'plot_range': (1,51,1), 
         'plot_levels':(0.1,0.3,1,3,10,30,100),
         'colors': ('white','lightgrey','grey','LightCoral','red','black')
         }

PCP_acc = {'name':'Precipitation (accumulated)',
           'shortname':'apcp', #'pcpacc',
           'grib_indicator':61, #gribtable['r'],
           'level_indicator': 105,
           'level': 457,#gribtable['acclevel'][0], 
           'tr_ind':0,#gribtable['intlevel'][1], # 0 or 4, 
           # test: new levels
           'grib_indicator_n':181, #gribtable['r'],
           'level_indicator_n': 105,
           'level_n': 0,#gribtable['acclevel'][0], 
           'tr_ind_n':4,#gribtable['intlevel'][1], # 0 or 4, 
           'grib_indicator_hir':61, #gribtable['r'],
           'level_indicator_hir': 105,
           'level_hir': 0,#gribtable['acclevel'][0], 
           'tr_ind_hir':4,
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30,100),
           'colors': ('white','lightgrey','grey','LightCoral','red','black'),
           #'plot_levels':(5,10,20,30,40,50,60,70,80,100),
           #'plot_levels':(20,40,60,80,100,120,140,160,180,200),
           #'colors': ('white','lightgrey','grey','LightCoral','red','lightblue','blue','darkblue','black'),
           # daily acc pictures on KNMI website:
           #'colors' : ((.992,.804,.027),(.796,.996,0.027),(.020,.804,.412),(.031,.702,.043),(.000,.525,.535),(.188,.302,.729)) 
           }

APCP_hir_tot = {'name':'Total precipitation (accumulated)',
           'shortname':'pcpacchir',
           'grib_indicator':61,
           'level_indicator': 105,
           'level': 0, 
           'plot_levels':(0.1,0.3,1,3,10,30,100),
           'colors': ('white','lightgrey','grey','LightCoral','red','black'),
           # daily acc pictures on KNMI website:
           #'colors' : ((.992,.804,.027),(.796,.996,0.027),(.020,.804,.412),(.031,.702,.043),(.000,.525,.535),(.188,.302,.729)) 
           }

SNOW_hir = {'name':'Snow',
           'shortname':'snow',
           'grib_indicator':79,
           'level_indicator': 105,
           'level': 0, 
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30),
           'colors': ('white','lightgrey','grey','LightCoral','red')
           }

NCPCP_hir = {'name':'Non-convective precipitation (accumulated)',
           'shortname':'ncpcphir',
           'grib_indicator':62,
           'level_indicator': 105,
           'level': 0, 
           'plot_levels':(0.1,0.3,1,3,10,30,100),
           'colors': ('white','lightgrey','grey','LightCoral','red','black'),
           # daily acc pictures on KNMI website:
           #'colors' : ((.992,.804,.027),(.796,.996,0.027),(.020,.804,.412),(.031,.702,.043),(.000,.525,.535),(.188,.302,.729)) 
           }

APCP_hir = {'name':'Convective precipitation (accumulated)',
           'shortname':'apcphir',
           'grib_indicator':63,
           'level_indicator': 105,
           'level': 0, 
           'plot_levels':(0.1,0.3,1,3,10,30,100),
           'colors': ('white','lightgrey','grey','LightCoral','red','black'),
           # daily acc pictures on KNMI website:
           #'colors' : ((.992,.804,.027),(.796,.996,0.027),(.020,.804,.412),(.031,.702,.043),(.000,.525,.535),(.188,.302,.729)) 
           }

SNOW    = {'name':'Snow',
           'shortname':'snow',
           'grib_indicator':62,#gribtable['s'],
           'level_indicator': 105,
           'level': 457,#gribtable['acclevel'][0], 
           'tr_ind':0,#gribtable['acclevel'][1], # 0 or 4 
           'grib_indicator_n':184,
           'level_indicator_n': 105,
           'level_n': 0,
           'tr_ind_n':4,
           'grib_indicator_hir':79,
           'level_indicator_hir': 105,
           'level_hir': 0,
           'tr_ind_hir':4,
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30),
           'colors': ('white','lightgrey','grey','LightCoral','red')
           }

SNOW_i    = {'name':'Snow',
           'shortname':'snow',
           'grib_indicator':62,#gribtable['s'],
           'level_indicator': 105,
           'level': 456,#gribtable['acclevel'][0], 
           'tr_ind':0,#gribtable['acclevel'][1], # 0 or 4 
           'grib_indicator_n':184,#gribtable['s'],
           'level_indicator_n': 105,
           'level_n': 0,#gribtable['acclevel'][0], 
           'tr_ind_n':0,# 0 or 4 
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30),
           'colors': ('white','lightgrey','grey','LightCoral','red')
           }


GRAUP   = {'name':'Graupel',
           'shortname':'graup',
           'grib_indicator':63,#gribtable['gr'],
           'level_indicator': 105,
           'level': 457, #gribtable['acclevel'], 
           'tr_ind':0, #gribtable['acclevel'][1], # 0 or 4, 
           'grib_indicator_n':201,#gribtable['s'],
           'level_indicator_n': 105,
           'level_n': 0,#gribtable['acclevel'][0], 
           'tr_ind_n':4,# 0 or 4 
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30),
           'colors': ('white','lightgrey','grey','LightCoral','red')
           }

GRAUP_i   = {'name':'Graupel',
           'shortname':'graup',
           'grib_indicator':63,#gribtable['gr'],
           'level_indicator': 105,
           'level': 456, #gribtable['acclevel'], 
           'tr_ind':0, #gribtable['acclevel'][1], # 0 or 4, 
           'grib_indicator_n':201,#gribtable['s'],
           'level_indicator_n': 105,
           'level_n': 0,#gribtable['acclevel'][0], 
           'tr_ind_n':0,# 0 or 4 
           #'plot_range': (1,51,1), 
           'plot_levels':(0.1,0.3,1,3,10,30),
           'colors': ('white','lightgrey','grey','LightCoral','red')
           }

SNWD    = {'name':'Snow depth',
           'shortname':'snowd',
           'grib_indicator':138,
           'level_indicator': 105,
           'level': 0,
           'plot_range': (1,51,1),
           'colormap': 'Blues'
            }
CLOUD   = {'name':'Cloud cover',
           'shortname':'tcc',
           'grib_indicator':71, #gribtable['tcc'],
           'level_indicator': 105,
           'level': 0, 
           'plot_range': (0,1.00001,1/8.), 
           'colormap':'binary_r' #'gray'
           }
MLCLOU  = {'name':'Cloud cover lowest model level',
           'shortname':'clwc',
           'grib_indicator':71,#gribtable['clwc'],
           'level_indicator': 109,
           'level': 60, 
           'plot_range': (0,1.1,1/8.), 
           'colormap':'gray'
           }
CLWAT   = {'name':'Cloud water',
           'shortname':'clwc',
           'grib_indicator': 76,#gribtable['clwc'],
           'level_indicator': 109,
           'level': 60, 
           'plot_range': (0,1.1,1/8.), 
           'colormap':'gray_r'
           }
VIS     = {'name':'Visibility',
           'shortname':'vis',
           'grib_indicator': 76,
           'level_indicator': 109,
           'level': 60, 
           'plot_levels': (0,20,50,100,200,400,700,1000,1500,2500,5e3,10e3,20e3), 
           'colors':('red','orange','yellow','GreenYellow','Lime','MediumSpringGreen',
                     'Cyan','RoyalBlue','MediumBlue','Indigo','DarkViolet','Fuchsia')
           }
CLICE   = {'name':'Cloud ice',
           'shortname':'ciwc',
           'grib_indicator': 58, #gribtable['ciwc'],
           'level_indicator': 109,
           'level': 60, 
           'plot_range': (0,1.1,1/8.), 
           'colormap':'gray_r'
           }
MLRA    = {'name':'Rain on lowest model level',
           'shortname':'pcpml',
           'grib_indicator': 62, #gribtable['r'],
           'level_indicator': 109,
           'level': 60, 
           'plot_levels': (0.1,0.3,1,3.,10,30), 
           'colormap':'Blues'
           }
MLSN    = {'name':'Snow on lowest model level',
           'shortname':'snowml',
           'grib_indicator': 79, #gribtable['s'],
           'level_indicator': 109,
           'level': 60, 
           'plot_levels': (0.1,0.3,1,3.,10,30), 
           'colormap':'Blues'
           }
MLGR    = {'name':'Graupel on lowest model level',
           'shortname':'graupml',
           'grib_indicator': 201, #gribtable['gr'],
           'level_indicator': 109,
           'level': 60, 
           'plot_levels': (0.1,0.3,1,3.,10,30), 
           'colormap':'Blues'
           }
SOILW   = {'name':'Soil water',
           'shortname':'soilw',
           'grib_indicator': 86, #gribtable['ssw'],
           'level_indicator': 105,
           'level': 801, 
           'plot_range': (0,0.5,0.05), 
           'draw_line':True,
           'line_range':(0,0.5,0.1),
           #'plot_levels': (0.1,0.3,1,3.,10,30), 
           'colormap':'Blues'
           }


G11 = T2M

#parameter_list = (T2M, RELH, CLOUD, U10, V10, PCP_i, PCP_acc, PRESS)
parameter_list = (T2M, CLOUD, PRESS,PCP_i, PCP_acc)

temp_list   = (T2M,RELH)
wind_list   = (U10, V10, UGST, VGST)

prec_list   = (PCP_i, PCP_acc, SNOW, GRAUP, CLWAT,
               MLRA, MLSN, MLGR, MLCLOU)#, 
               #APCP_hir,NCPCP_hir,MLRA, MLSN, MLGR, MLCLOU)

prec_list_hir = (APCP_hir_tot,SNOW_hir,MLRA,MLCLOU)

cloud_list  = (CLOUD,BRTEMP)

point_list  = (PCP_acc,T2M,)

postpr_list = (T2M, RELH, 
               PCP_i, SNOW, GRAUP, MLRA, MLSN, MLGR
               )

all_list = (default,T2M, RELH, PRESS, 
            CLOUD, CLWAT, CLICE,  
            U10, V10, UGST, VGST, WIND, VV,
            PCP_i, PCP_acc, GRAUP, SNOW, SNWD, MLRA, MLSN, MLGR,
            APCP_hir_tot,SNOW_hir,
            SOILW,VIS)

#return parameter_list
