#! /usr/bin/env python

'''
A command line interface to the PARSE_VOBSVFLD.py module

usage:


'''

import os,sys,argparse

import parse_vobsvfld as pv

## 
#if __name__ == '__main__':

parser = argparse.ArgumentParser(description='parse vobs/vfld files to SQLite database (in outdir), inputfiles of type vfld or vobs allow posix wildcards, version is detected')
parser.add_argument('-V','--version',
                    action='version',
                    version='0.1')
parser.add_argument('-v', '--verbose',
                    default=False,
                    action='store_true',
                    help="Give more verbose output." )
parser.add_argument('-o', '--outdir',
                    action='store',
                    dest='outdir', 
                    required=True,
                    help="The ouput directory" )
parser.add_argument('infiles', metavar='inputfiles', type=str,      action='append',     nargs='*',
                    help='files to parse')

args = parser.parse_args()

# check input file
if not args.infiles:
    assert False, "no input files"

flattend = []
if args.infiles:
    for sublist in args.infiles:
        flattend += sublist 
    args.infiles = flattend

# check outputfile
if not args.outdir:
    assert False, "no output directory given"    

if args.verbose:
    print 'Input files:'
    for f in args.infiles: print '>> ',f
    print 'Output to ',args.outdir

# call process_v(vfile,outdir = './',tempdir=None,model=None,vdate=None)
# or parse_v(vobsfile,vdate=None,model='default',leadtime=0,outdir = './')
for f in args.infiles:
    pv.process_v(f,outdir = args.outdir)
