#! /usr/bin/env python

import os,sys
import datetime as dt
import subprocess

def get_remote(subdir,outdir,rfile = None,rexp = 'ex*',obrexp = 'vfld*gz'):
    '''
    Copies output of extract_remote() (above) and obs files to local computer 
    to perform verification (etc) 
    '''
    
    print 'At ',dt.datetime.today()
    dsa = '/usr/people/plas/.ssh/id_dsa'

    if rfile == None:
        ccmd = 'scp -i '+dsa+' testharm@bxshpc01.knmi.nl:/nfs_ltc/lam_test/testharm/pcptemp/'+subdir+'/'+rexp+' '+outdir
    else:
        ccmd = 'scp -i {dsa} testharm@bxshpc01.knmi.nl:/nfs_ltc/lam_test/testharm/pcptemp/{sub}/{rf} {local}'.format(dsa=dsa,sub=subdir,rf=rfile,local=outdir)
    print ccmd
    subcmd = subprocess.Popen(ccmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    #check: stdout = subprocess.check_output(ccmd+"; exit 0",stderr=subprocess.STDOUT,shell=True)
    # instead of os.system(ccmd)

    ocmd = 'scp -i '+dsa+' testharm@bxshpc01.knmi.nl:/nfs_ltc/lam_test/testharm/pcptemp/'+subdir+'/'+obrexp+' '+outdir
    print ocmd
    #os.system(ocmd)

    return subcmd.stderr.read() #stdout

