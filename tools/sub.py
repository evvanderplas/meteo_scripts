#! /usr/bin/env python

import os,sys

if __name__ == '__main__':

    args = sys.argv
    print 'in sub: ',args, os.environ['SUBSTRING']

    if len(args) > 1:
        outf = os.path.join('subout','file_{nr}.txt'.format(nr = args[1]))
        print outf
        f = open(outf,'w')
        f.write('Just wrote something in file nr {nr}'.format(nr = args[1]))
        f.close()

        # some failing code: 
        fail = open('nonexisting.txt','r')
        try:
            fail = open('nonexisting.txt','r')
        except:
            print 'something failed'
