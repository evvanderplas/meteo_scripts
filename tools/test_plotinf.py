#! /usr/bin/env python

import os,sys

import numpy as np
import datetime as dt
import pygrib
#import h5py

class plotinfo():

    def __init__(self,model='Harmonie',arg1=None,arg2=None,arg3='drie',**kwargs):
        
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3

        print kwargs
        for k,v in kwargs.items():
            print k,v
            self.__dict__[k] = v
    
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
    
    my_inst = plotinfo(**kwargs)

    # test filetype, assume it is grib if there is no extension:
    if os.path.splitext(datafile)[1] in ('.h5','.hdf5'):
        print 'hdf5' #my_inst.read_from_h5(datafile) #,**kwargs)
    elif os.path.splitext(datafile)[1] in ('.grb','.grib',''):
        print 'grib' #my_inst.read_from_grib(datafile) #,**kwargs)

    return my_inst

if __name__ == '__main__':

    print 'testje'
    mypl = read('/not/a/file.h5',arg3='vier',onzin='ja',nounou='soms',getal = 3)
