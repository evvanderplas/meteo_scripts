#!/usr/bin/env python

"""
a litle module that allows inspecting size and growth of memory usage of the
current python process. Usefull for tracing memory leaks.

Written by: J. de Kloe, 2012
"""

import os #,sys
import numpy

prev_vsize      = None

def report_mem_usage():
    """
    method: read stat info from the /proc filesystem and extract
    the appropriate field (works for linux only)
    """
    global prev_vsize

    statfile = '/proc/self/stat'
    txt = open(statfile).read()
    parts = txt.split()
    #print parts
    #print 'len = ',len(parts)
    # for field definitions see: man 5 proc
    # or: http://www.kernel.org/doc/man-pages/online/pages/man5/proc.5.html
    #
    # interesting for this script:
    # 22: vsize : virtual memory size in bytes

    vsize      = int(parts[22])

    if prev_vsize:
        print 'vsize change:      ', vsize - prev_vsize

    prev_vsize      = vsize

if __name__ == "__main__":
    # first call is to init prev_size field
    report_mem_usage()
    
    print '--- 50000 zeros ---'
    x = numpy.zeros(50000)
    report_mem_usage()

    print '--- 100000 zeros ---'
    y = numpy.zeros(100000)
    report_mem_usage()

    print '--- 200000 zeros ---'
    z = numpy.zeros(200000)
    report_mem_usage()

    print '--- single character ---'
    a1='a'
    report_mem_usage()

    print '--- two characters ---'
    a2='aa'
    report_mem_usage()

    print '--- done ---'
