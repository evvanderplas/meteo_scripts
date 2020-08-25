#! /usr/bin/env python

'''
Script that tests the generation of lats,lons from MSG satellite images. 
One is a SAFNWC MSGCPP product, which is cut a liitle differently from the 
raw image that is retrieved from EUMETSAT. 

A projection string is fed to pyproj, and with the corners as calculated from the MSG 
metadata we can get lats,lons from this projection object using its inverse method. 

'''

import os,sys,glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

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

td = dt.datetime.today()
td0 = td.replace(hour = td.hour -3,minute=0,second=0,microsecond=0)
dtg= {'y':td0.year,
      'm':str(td0.month).zfill(2),
      'd':str(td.day).zfill(2),
      'h':str(td.hour).zfill(2),
      }
#print td0; sys.exit(1)


grbirfile = '/nobackup_1/users/plas/temp/HARM_N55_201210250000_01200_GB'
msgirfile = os.path.join(MSGIR['dir'],'METEOSAT_10_SEVIRI_EUROPE_2014_02_03_07_00_00.h5')
msgsaffile = os.path.join(MSG['dir'].format(y=td0.year,m= str(td0.month).zfill(2),d=str(td.day).zfill(2) ),MSG['fformat'].format(dtgh=td0.strftime('%Y%m%d%H')))
print grbirfile, os.path.exists(grbirfile)
print msgirfile, os.path.exists(msgirfile)
print msgsaffile, os.path.exists(msgsaffile)
#sys.exit(0)
            
## msg METEOSAT_10 proj:
'''
geo_column_offset = -1500.0
geo_row_offset = -1800.0
geo_pixel_size_x,geo_pixel_size_y = 3.000387,3.000387

x = (3000 + geo_column_offset) * geo_pixel_size_x
y = (1000 + geo_row_offset) * geo_pixel_size_y

Maarten Plieger:
gdal_translate -of G  -a_ullr 0 -3649.998 700.0007 -4415.002 -a_srs "+proj=stere +lat_0=90 +lon_0=0 +lat_ts=60 +a=6378.14 +b=6356.75 +x_0=0 y_0=0" in.h5 out.tif 
gdal_translate -of G  -a_ullr x_ul y_ul x_lr y_lr -a_srs sourceprojstring in.h5 out.tif 

'''

x_lr =   4500.5805
x_ul =  -4500.5805 
y_lr =   2400.3096
y_ul =   5400.6966

projs = '+proj=geos +a=6378.1690 +b=6356.5838 +lat_0=0.0 =+lon_0=0.0 +h=35785.8310'

xs = np.linspace(x_ul, x_lr, 3001)
ys = np.linspace(y_ul, y_lr, 1001)
xn,yn = np.meshgrid(xs,ys)

import pyproj
p = pyproj.Proj(projs)
lons,lats = p(xn,yn,inverse=True)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.contourf(lons,np.arange(-90,90,5))
ax2.contourf(lats,np.arange(0,90,5))
ax1.set_title('Longitudes')
ax2.set_title('Latitudes')

fname1 = 'MSG_proj_1.png'
fig.savefig(fname1)
print 'Created ',fname1

from mpl_toolkits.basemap import Basemap # map functionality
maph = Basemap(projection = 'geos',
               rsphere=(6378.16900,6356.5838),
               lon_0=0.0,lat_0=0.0,
               llcrnrx = x_ul,llcrnry = y_lr,
               urcrnrx = x_lr,urcrnry = y_ul,
               satellite_height=35785.8310
               )

mlons, mlats  = maph(xn,yn,inverse=True)

saf = True
if saf: # from hdf5 file
    x_lr =  5565.7482275008615
    x_ul = -5568.748630858005
    y_lr =  2655.3569710718393
    y_ul =  5436.730883143699

    xs = np.linspace(x_ul, x_lr, 3712)
    ys = np.linspace(y_ul, y_lr, 928)
    xn,yn = np.meshgrid(xs,ys)
    slons,slats = p(xn,yn,inverse=True)

    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax1.contourf(slons,np.arange(-90,90,5))
    ax2.contourf(slats,np.arange(0,90,5))
    ax1.set_title('Longitudes')
    ax2.set_title('Latitudes')
    
    fname2 = 'MSG_proj_2.png'
    fig.savefig(fname2)
    print 'Created ',fname2
