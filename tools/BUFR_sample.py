#! /usr/bin/env python

import os,sys

import pybufr_ecmwf
#  #[ imported modules
#import sys # operating system functions
#import getopt # a simpler version of argparse, which was introduced in
              # python 2.7, and is not by default available for older versions

# import the python file defining the RawBUFRFile class
from pybufr_ecmwf.bufr import BUFRReader
#from pybufr_ecmwf.raw_bufr_file import RawBUFRFile
#from pybufr_ecmwf.bufr_interface_ecmwf import BUFRInterfaceECMWF
#from pybufr_ecmwf.helpers import python3

#  #]


def print_bufr_content1(input_bufr_file, output_fd, separator, max_msg_nr):
    #  #[ implementation 1
    """
    example implementation using the BUFRReader class
    combined with the get_values_as_2d_array method
    """
    
    # get an instance of the BUFR class
    # which automatically opens the file for reading and decodes it
    bob = BUFRReader(input_bufr_file, warn_about_bufr_size=False)

    msg_nr = 0
    while True:
        try:
            bob.get_next_msg()
            msg_nr += 1
        except EOFError:
            break

        # add header strings
        # print 'DEBUG: bob.msg_loaded ',bob.msg_loaded
        if bob.msg_loaded == 1:
            list_of_names = []
            list_of_units = []
            list_of_names.extend(bob.get_names())
            list_of_units.extend(bob.get_units())
            # print 'DEBUG: ',separator.join(list_of_names)
            # print 'DEBUG: ',separator.join(list_of_units)
            output_fd.write(separator.join(list_of_names) + "\n")
            output_fd.write(separator.join(list_of_units) + "\n")
        
        data = bob.get_values_as_2d_array()
        # print 'DEBUG: data.shape = ', data.shape
        if data.shape[0]*data.shape[1] == 0:
            print 'NO DATA FOUND! this seems an empty BUFR message !'
            continue

        for subs in range(len(data[:, 0])):
            output_fd.write(str(subs)+separator+
                            separator.join(str(val) for val in data[subs, :])+
                            "\n")
        print 'converted BUFR msg nr. ', msg_nr
        if ( (max_msg_nr>0) and (msg_nr >= max_msg_nr) ):
            print 'skipping remainder of this BUFR file'
            break

    # close the file
    bob.close()
    if msg_nr == 0:
        print 'no BUFR messages found, are you sure this is a BUFR file?'
    #  #]


if __name__ == '__main__':

    bufrdir = '/net/bhw379/nobackup/users/plas/temp/BUFRdata/'
    bufrfile = os.path.join(bufrdir,'BUFR.radarv.bewid')
    os.environ['BUFR_TABLES'] = os.path.join(bufrdir,'BUFRtables')


    bob = BUFRReader(bufrfile, warn_about_bufr_size=False)
    bob.get_next_msg()
    names = bob.get_names()
    units = bob.get_units()
    
    print set(names)
    print 20*'='
    print set(units)
    for n,u in zip(names,units):
        if "dB" in u: print n,u


    outf = open('radarBUFR.dat','w')
    print_bufr_content1(bufrfile, outf, ',' , 5)
