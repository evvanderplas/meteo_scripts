#! /usr/bin/env python

import os,sys,plotinfo

msgfile = '/data/mos/obsmsg/seviri/level2/safnwc/2013/10/15/SAFNWC_MSG3_CMa__201310151800_MSG-N_______.h5'
msgfile = 'SAF_MSG3_CMa_201308060400.h5'
print msgfile, os.path.exists(msgfile)
m = plotinfo.plotinfo(model='msg')
m.read_from_h5(msgfile,datapath='CMa',name='MSG')
ll = m.latlons
dim0,dim1 = ll[0].shape
k,l = (dim0/2) - 300, (dim1/2) - 300
m,n = (dim0/2) + 300, (dim1/2) + 300
print 'median lat,lon:',ll[0][dim0/2,dim1/2],ll[1][dim0/2,dim1/2]

ulats,ulons = ll[0][k:m,l:n],ll[1][k:m,l:n]
print 'SE,NW:',ulats[0,0],ulons[0,0],' ',ulats[-1,-1],ulons[-1,-1]
sys.exit(0)

home = os.environ['HOME']
datadir  = os.path.join(home,'HARP/Harp_sample/data')
gribfile = os.path.join(datadir,'fc','harm','ex2012062100+012')
gribfile = 'ex2013080600+004grib_sa'
gr12 = plotinfo.plotinfo(model='harm')
gr12.read_from_grib(gribfile,71,0,model = 'test')

m.resample(gr12)
m.plot(domain='nl')
