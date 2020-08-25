#! /usr/bin/env python

'''
sample script to demonstrate use of bmap functionality
'''

import os,sys,glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

levels,colors = None,None

MSGIR  = {'model':'MSG',
          'name':'MSG IR',
          'shortname':'MSGIR',
          'dir':'/net/bhw422/nobackup/users/valkde/MSG-hdf', # /yyyy/mm/dd/
          'rexp':'METEOSAT_10_SEVIRI_EUROPE_*00_00.h5',
          'fformat':'METEOSAT_10_SEVIRI_EUROPE_{y}_{m}_{d}_{h}_00_00.h5',
          'datapath':'image5/image_data',
          'metapath':'image5/satellite',           
          }
MSG    = {'model':'MSG',
          'name':'MSG Cloud mask (CMa)',
          'shortname':'CM',
          'dir':'/data/mos/obsmsg/seviri/level2/safnwc/{y}/{m}/{d}/',
          'rexp':'SAFNWC_MSG2_CMa_*00_MSG*.h5',
          'fformat':'SAFNWC_MSG3_CMa__{dtgh}00_MSG-N_______.h5', #SAFNWC_MSG3_CMa__201401311000_MSG-N_______.h5
          'datapath':'CMa',
          'metapath':'',
          }

radfile = '/nobackup/users/plas/radar/RAD_NL25_PCP_NA_201008131100.h5' 
radfile = './RAD_NL25_PCP_NA_201312051600.h5' 
print radfile, os.path.exists(radfile); #sys.exit(0)

gfile = '/nobackup/users/plas/metview/exp/front/36h12/20100813/hirp0/fc2010081300+0720grib'
gfile     = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013123000+015grib_sa'
gfileprev = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013123000+014grib_sa'

gfile     = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013061900+012grib'
#gfileprev = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013061900+007grib_sa'

gfile     = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013120500+016grib_sa'
gfileprev = '/nobackup_1/users/plas/verif/BULL/BULL/ex2013120500+015grib_sa'

td = dt.datetime.today()
td0 = td.replace(hour = td.hour -3,minute=0,second=0,microsecond=0)
dtg= {'y':td0.year,
      'm':str(td0.month).zfill(2),
      'd':str(td.day).zfill(2),
      'h':str(td.hour).zfill(2),
      }
#print td0; sys.exit(1)


grbirfile = '/nobackup_1/users/plas/temp/HARM_N55_201210250000_01200_GB'
msgirfile = os.path.join(MSGIR['dir'],MSGIR['fformat'].format(y=td0.year,m= str(td0.month).zfill(2),d=str(td.day).zfill(2),h= str(td0.hour).zfill(2) ))
msgsaffile = os.path.join(MSG['dir'].format(y=td0.year,m= str(td0.month).zfill(2),d=str(td.day).zfill(2) ),MSG['fformat'].format(dtgh=td0.strftime('%Y%m%d%H')))
print grbirfile, os.path.exists(grbirfile)
print msgirfile, os.path.exists(msgirfile)
print msgsaffile, os.path.exists(msgsaffile)


print gfile, os.path.exists(gfile)

import bmap
#import sources, settings
import plotinfo as getinfo

gr12 = getinfo.plotinfo(model='model')
gr12.read_from_grib('samplegrib_LA.grb',61,level=457,TR=0,case='penalty')

if 1:
    print 'opening ',msgirfile,os.path.exists(msgirfile)
    msgir07 = getinfo.plotinfo(model='msgir')
    msgir07.read_from_h5(msgirfile,datapath=MSGIR['datapath'],metapath=MSGIR['metapath'],model=MSGIR['shortname'],name=MSGIR['name']) #msgir07.get_settings(preset = '')

    #sys.exit(0)

    msgir07.resample(gr12,nodata=65535) 
    msgir07.plot_levels = np.linspace(-40,-10,31)
    msgir07.plot(domain= 'eur',lsmask_colour = 'red')  

    msg = getinfo.plotinfo(model='msg')
    msg.read_from_h5(msgsaffile,datapath=MSG['datapath'],metapath=MSG['metapath'],model=MSG['shortname'],name=MSG['name'])
    msg.resample(gr12,nodata=65535) 
    msg.plot_levels = np.linspace(0,1,11)
    msg.plot(domain= 'eur',lsmask_colour = 'red')  

    gir = getinfo.plotinfo(model='harmonie')
    gir.read_from_grib(grbirfile,118,level=39680,TR=0,case='penalty')

    sys.exit(0)

if 1:
    print 'opening ',radfile,os.path.exists(radfile)
    rad12 = getinfo.plotinfo(model='radar')
    rad12.read_from_h5(radfile)
    rad12.get_settings(preset = 'APCP')
    #rad12.plot(domain = 'nl')
    
if 1:
    print 'opening ',gfile,os.path.exists(gfile)
    gr12 = getinfo.plotinfo(model='model')
    if 0:
        gr12.read_from_grib(gfile,61,level=456,TR=0,case='penalty')
        gr12.values = 3600 * gr12.values
    else:
        gr12.read_from_grib(gfile,61,level=457,TR=0,case='penalty')

        gr11 = getinfo.plotinfo(model='model')
        gr11.read_from_grib(gfileprev,61,level=457,TR=0,case='penalty')

        gr12.values = gr12.values - gr11.values

    rad12.resample(gr12)
    gr12.get_settings(preset = 'PCP')
    #gr12.plot(domain = 'nl')

#sys.exit(0)
fig  = plt.figure()
ax   = fig.add_subplot(1,1,1)
mplot = bmap.myMap(domain='nl',dbg = True)
mplot.xymap(latlons = gr12.latlons, domain = 'nl')
mplot.set_axes(ax)
mplot.dress_up_map(domain = 'nl',lsmask_colour = 'k')
gcont = mplot.bmap.contourf(mplot.x,mplot.y,gr12.values,
                            [2.,100],colors = 'red',ax=ax)
rcont = mplot.bmap.contourf(mplot.x,mplot.y,rad12.values,
                            [2.,100],colors = 'darkgreen',ax=ax)

ax.set_title('Forecast (red) and radar (green) \n for heavy rain (>2mm/h), 5 Dec 2013, 16:00 UTC')
fig.savefig('custom.png')
