#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle
import time
import datetime as dt
import numpy
import matplotlib, matplotlib.pyplot
from mpl_toolkits.basemap import Basemap # map functionality
from matplotlib import cm as CM
import pygrib

fig  = matplotlib.pyplot.figure()
ax   = fig.add_subplot(1,1,1)

f = open('nlmap.pkl','r')
bmap = pickle.load(f)
f.close()
bmap.drawcoastlines(ax=ax)
bmap.drawparallels([52,53])
mrd = bmap.drawmeridians([5,6])

grbs = pygrib.open('./brtemp.grb')
w          = grbs[1].values
lats, lons = grbs[1].latlons()
x,y        = bmap(lons,lats)


jcmap   = CM.get_cmap('Reds')
#cmap.set_over( color = 'black') #,alpha = 0.5)
jcmap.set_over( 'black')
jcmap.set_under('black')

print jcmap(-0.3)
print jcmap(0.)
print jcmap(1.)
print jcmap(1.1)

levels = numpy.arange(-255,-220,5)

#sys.exit(1)
print w.min(), w.max()

g  = bmap.contourf(x,y,w,levels,cmap = jcmap,extend='max',ax = ax)
cb = fig.colorbar(g,shrink=0.9,format='%.1f') # '%.1e'

print g.cmap(-0.3)
print g.cmap(0.)
print g.cmap(1.)
print g.cmap(1.1)

bf =  g.cmap._lut
g.cmap.set_over( 'black')
g.cmap.set_under(color = 'k')
#cb.Normalize(clip=False)
aft =  g.cmap._lut
print numpy.all(bf == aft )

sf = fig.savefig('bmap1.png')
print 'created plot: bmap1.png'

#print mrd
#print dir(g.cmap)
#for item in dir(g):
#    print item, type(getattr(g,item))
#print g.tcolors

del(mrd[5])
sf = fig.savefig('bmap2.png')
print 'created plot: bmap2.png'


#print dir(bmap)
