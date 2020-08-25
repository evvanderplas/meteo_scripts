#!/usr/bin/env python

"""
a little script to install the pygrib python library and its
dependencies on fedora
"""

#  #[ method
# Download latest versions of the needed source codes (replace the crosses
# by the most current version numbers):
#   * grib_api-x.x.x.tar.gz from: http://www.ecmwf.int/products/data/\
#                                        software/download/grib_api.html
#   * pyproj-x.x.x.tar.gz from: http://code.google.com/p/pyproj/downloads/list
#   * pygrib-x.x.x.tar.gz from: http://code.google.com/p/pygrib/downloads/list
#
# huidige versies:
#   grib_api-1.9.9.tar.gz
#   pyproj-1.8.9.tar.gz
#   pygrib-1.8.3.tar.gz
#
# volg dan de instructies op:
# http://code.google.com/p/pygrib/wiki/LinuxMacInstallation
#
#
# * mkdir grib_api_dir
# * tar zxvf grib_api-x.x.x.tar.gz
# * cd grib_api-x.x.x
# * (for bash) export CFLAGS="-O2 -fPIC" 
# * (for c-shell) setenv CFLAGS "-O2 -fPIC"
# * ./configure --prefix=`pwd`/../grib_api_dir
# * make
# * make check
# * make install
# * cd ..
#
# * tar zxvf pyproj-x.x.x.tar.gz
# * cd pyproj-x.x.x
# * python setup.py install --user
# * cd ..
#
# * (for bash) export GRIBAPI_DIR=`pwd`/grib_api_dir
# * (for bash) export JASPER_DIR=/usr
# * (for c-shell) setenv GRIBAPI_DIR `pwd`/grib_api_dir
# * (for c-shell) setenv JASPER_DIR /usr
# * tar zxvf pygrib-x.x.x.tar.gz
# * cd pygrib-x.x.x
# * python setup.py install --user
# * python test.py
#  #]
#  #[ imported modules
import os, sys, glob
from fabric.api import local, lcd, settings
# for fabric docs see: http://docs.fabfile.org/en/1.0.1/index.html
#  #]
#  #[ settings
install_grib_api = True
install_pyproj   = True
install_pygrib   = True

source_dir         = '/usr/local/free/src/'

grib_api_version = '1.11.0'
version_pyproj   = '1.9.3'
version_pygrib   = '1.9.7'
fedora_version   = 'f18'
# don't use this one! gives this error
#     import pygrib
#     ImportError: /usr/local/free/lib/python2.7/site-packages/pygrib.so:
#     undefined symbol: opj_setup_decoder
#version_pygrib   = '1.8.4'

set_permissions_tool = '/net/bhw412/nobackup/users/kloedej/'+\
                       'versiebeheer_mercurial/bin_scripts/file_tools/'+\
                       'set_permissions_insfree_rx.py'

install_dir_grib_api = '/usr/local/free/installed/grib_api-'+\
                       grib_api_version+'-PIC-for-pygrib-'+fedora_version
build_dir_grib_api = source_dir+'pygrib-install'
grib_tarfile = 'grib_api-'+grib_api_version+'.tar.gz'
grib_source_dir = source_dir+'grib_api-'+grib_api_version+'-src-for-pygrib'

pyproj_tarfile = 'pyproj-'+version_pyproj+'.tar.gz'
pyproj_source_dir = source_dir+'/pyproj-'+version_pyproj

pygrib_tarfile = 'pygrib-'+version_pygrib+'.tar.gz'
pygrib_source_dir = source_dir+'/pygrib-'+version_pygrib

# installation in: /usr/local/free/lib/python2.7/site-packages/
install_dir_pyproj = '/usr/local/free'
install_dir_pygrib = '/usr/local/free'
#  #]
if install_grib_api:
    #  #[
    if os.path.exists(grib_source_dir):
        print 'Deleting old source dir: ',grib_source_dir
        os.system('\\rm -rf '+grib_source_dir)

    local('mkdir '+grib_source_dir)
    with lcd(grib_source_dir):
        # this check is to see if the cd command really works
        # Note that Fabric-1.0.1-py2.7.egg which I had installed locally for
        # some time has a broken cd command and cannot be used
        # with this script!!!
        local('pwd')
        local('tar zxvf ../'+grib_tarfile)

    # werken met env settings gedefinieerd met het 'with' statement
    # werkt nog niet goed met de hudige fabric versie:
    # with settings(cd(grib_source_dir+'/grib_api-'+grib_api_version),
    #               CFLAGS="-O2 -fPIC"):
    #        local('./configure --prefix='+install_dir_grib_api)
    # er was al ee bugrapport hiervoor.
    # zie: http://code.fabfile.org/issues/show/316
    #  en: http://code.fabfile.org/issues/show/252
    #
    with lcd(grib_source_dir+'/grib_api-'+grib_api_version):
        local('export CFLAGS="-O2 -fPIC";'+\
              './configure --enable-python --prefix='+install_dir_grib_api)

    with lcd(grib_source_dir+'/grib_api-'+grib_api_version):
        local('make')

        # note: the make check gives an error for version 1.9.9, therefore
        # this warn_only setting is required!
        with settings(warn_only=True):
            r=local('make check')
            if r.failed: print r
            r=local('make install')
            if r.failed: print r

    cmd = set_permissions_tool+' '+grib_source_dir
    os.system(cmd)
    cmd = set_permissions_tool+' '+install_dir_grib_api
    os.system(cmd)    
    #  #]

if install_pyproj:
    #  #[
    if os.path.exists(pyproj_source_dir):
        print 'Deleting old source dir: ',pyproj_source_dir
        os.system('\\rm -rf '+pyproj_source_dir)

    with lcd(source_dir):
        local('tar zxvf '+pyproj_tarfile)
    with lcd(pyproj_source_dir):
        # local('python setup.py install --user')
        # of
        local('python setup.py install --prefix=/usr/local/free')

    cmd = set_permissions_tool+' '+pyproj_source_dir
    os.system(cmd)
    cmd = set_permissions_tool+' '+\
          glob.glob('/usr/local/free/lib*/python2.7/site-packages/pyproj')[0]
    os.system(cmd)
    cmd = set_permissions_tool+' '+\
          glob.glob('/usr/local/free/lib*/python*/site-packages/'+\
                    'pyproj-*.egg-info')[0]
    os.system(cmd)
    #  #]

if install_pygrib:
    #  #[
    if os.path.exists(pygrib_source_dir):
        print 'Deleting old source dir: ',pygrib_source_dir
        os.system('\\rm -rf '+pygrib_source_dir)

    with lcd(source_dir):
        local('tar zxvf '+pygrib_tarfile)
    
    with lcd(pygrib_source_dir):
        env1 = 'setenv GRIBAPI_DIR '+install_dir_grib_api    
        env2 = 'setenv JASPER_DIR /usr'
        env3 = 'setenv OPENJPEG_DIR /usr'
        env4 = 'setenv PNG_DIR /usr'
        env5 = 'setenv ZLIB_DIR /usr'
        cmd = 'python setup.py install --prefix=/usr/local/free'
        local('csh -c "'+env1+';'+env2+';'+env3+';'+env4+';'+env5+';'+\
              cmd+'"')

    cmd = set_permissions_tool+' '+pygrib_source_dir
    os.system(cmd)

    #cmd = 'python test.py'
    #os.system(cmd)
    
    cmd = set_permissions_tool+' '+\
          glob.glob('/usr/local/free/lib*/python*/site-packages/pygrib.so')[0]
    os.system(cmd)
    cmd = set_permissions_tool+' '+\
          glob.glob('/usr/local/free/lib*/python*/site-packages/g2clib.so')[0]
    os.system(cmd)
    cmd = set_permissions_tool+' '+\
          glob.glob('/usr/local/free/lib*/python*/site-packages/'+\
                    'pygrib-*.egg-info')[0]
    os.system(cmd)

# run unit tests for pygrib:
#* python test.py
#
# note, for now this env setting is needed !
#setenv PYTHONPATH /usr/local/free/lib/python2.7/site-packages/
    #  #]
