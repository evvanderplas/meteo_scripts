#! /usr/bin/env python

from harmonie_codes import *

w = 1
if w:
    f = open('harmonie_codes36.txt','w')
    f.write('field_codes_harmonie36 = {\n')

for k in sorted(field_codes_harmonie37.keys()):
    nk = "'"+'001'+k[3:]+"'"
    v  = field_codes_harmonie37[k]
    #print ''.join(['    ',nk,':',str(field_codes_harmonie37[k]),','])
    ln,sn,u,fa = "'"+v[0]+"'","'"+v[1]+"'","'"+v[2]+"'","'"+v[3]+"'"
    try:
        stn = "'"+v[4]+"'"
    except:
        stn = "'None'"
    hcode = "    {nk} : [{ln:<35},{sn:<15}, {u:<15},{fa:<20},{stn:<30}],".format(nk=nk,ln=ln,sn=sn,u=u,fa=fa,stn=stn)
    print hcode
    if w: f.write(hcode+'\n')
if w:
    f.write('}')
    f.close() 
