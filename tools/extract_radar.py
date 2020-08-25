#! /usr/bin/env python

import os,sys,glob
import datetime as dt

def suf2(h):
    '''take int,float or string, add 0 if smaller than 10, return string: 1 -> '01' etc '''
    if int(h) < 10: return '0'+str(h)
    else: return str(h)

def locate(pattern, root=os.curdir):
    '''
    Locate all files matching supplied filename pattern in and below
    supplied root directory.
    '''
    import fnmatch

    for path, dirs, files in os.walk(os.path.abspath(root)):
        #print files
        for filename in fnmatch.filter(files, pattern):
            #print filename
            yield os.path.join(path, filename)

def extract_radar(radar_archive,outdir,tempdir=None,rexp = 'RAD_NL23*',debug=True):

    '''Extracting only relevant (eg RAD_NL23*) files from tar tape archive '''

    import tarfile,shutil

    def rrad_files(members):
        '''adapted from python website'''
        if '*' in rexp.strip('*'):
            r = rexp[:rexp.find('*')]
        else:
            r = rexp.strip('*')
        for tarinfo in members:
            #mrexp = rexp[:rexp.find('*')] #;sys.exit(1)
            #print rexp,rexp.strip('*'); sys.exit(1)
            if r in os.path.split(tarinfo.name)[1]:
                yield tarinfo

    filename = os.path.split(radar_archive)[1]
    print '... extracting radar files ...', filename,'with pattern',rexp

    ar = tarfile.open(radar_archive,'r:*') # open for reading with transparent compression
    #m = ar.getmembers()

    if tempdir == None:
        tempdir = os.path.join(outdir,'temp')
    elif os.path.isdir(tempdir):
        pass
    elif not os.path.isdir(tempdir):
        try:
            os.makedirs(tempdir)
        except OSError:
            print 'Dir exists?',tempdir
            #print 'Hm, no dir, not makeable: ',tempdir; sys.exit(1)

    print 'untarring to ',tempdir
    outsafe = tempdir #os.path.join(outdir,'temp') # os.path.join(outdir,suf2(dtf.year),suf2(dtf.month))
    print filename,' to: ',outsafe

    if debug:
        print 'what in tarfile?'
        ar.list()
        for m in rrad_files(ar): print m

    l = ar.extractall(members=rrad_files(ar),path=outsafe)
    print 'using ', rexp,' in locate()',outsafe 
    for xfile in locate(rexp, root=outsafe):
        print xfile
        try:
            shutil.move(xfile,outdir)
        except:
            os.rename(xfile,os.path.join(outdir,os.path.split(xfile)[1]))

    #sys.exit(1)
    return 'done'

if __name__ == '__main__':


    dtg  = dt.datetime(2012,3,7,0)
    dte  = dt.datetime(2012,3,7,23)

    dtg  = dt.datetime(2013,12,5,0)


    tarrexp  = 'OMBE_HDF5_RADAR_*20120307*.tar'
    
    ## radar archive directory
    tardir   = os.path.join('/data/mos/fa/ao/ombe/ombe_lt/radar/archive/hdf',str(dtg.year))
    tardir   = os.path.join('/data/mos/obsombe/ombe_lt/radar/archive/hdf',str(dtg.year))
    tarfiles = os.path.join(tardir,'OMBE_HDF5_RADAR_{dtg}.tar'.format(dtg=dtg.strftime('%Y%m%d')))

    for sfile in [tarfiles,]: 
        print os.path.split(sfile)[1]
        extract_radar(sfile,'./',rexp = 'RAD_NL25*')

    sys.exit(1)
