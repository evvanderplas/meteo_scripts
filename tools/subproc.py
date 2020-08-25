#! /usr/bin/env python

##
## script to test some subprocess functionality
## 

import os,sys,subprocess

if __name__ == '__main__':

    
    for i in range(10):
        print 'starting process ',i
        os.environ['SUBSTRING'] = 'sub'+str(i)
        process = subprocess.Popen(['python', 'sub.py', str(i)], stdout=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print 'stdout: ',stdout
        print 'stderr: ',stderr
        
