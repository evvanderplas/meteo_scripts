#! /usr/bin/env python

import os,sys,sqlite3
from os.path import join as pjoin

dbdir = '/nobackup/users/plas/harp'
fformat = 'synop_{m}_{p}_{ym}.db'

m = 'RUC'
#m = 'H4DVAR'
#m = 'BULL'
subdir = 'ecjan'
#subdir = 'BULL_3H'
p = 'wind'
ym = '201402'
db = pjoin(dbdir,subdir,fformat.format(m=m,p=p,ym=ym))

print db, os.path.exists(db)

if os.path.exists(db):
    conn   = sqlite3.connect(db)
    c      = conn.cursor()

    sql    = 'select * from {p} where sid between 6000 and 7000'.format(p=p)
    c.execute(sql)
    data = c.fetchall()
    conn.close()

p2 = 'winddir'
db2 = pjoin(dbdir,subdir,fformat.format(m=m,p=p2,ym=ym))
print db2, os.path.exists(db2)

if os.path.exists(db2):
    conn2   = sqlite3.connect(db2)
    c2      = conn2.cursor()

    sql    = 'select * from {p} where sid between 6000 and 7000'.format(p=p2)
    c2.execute(sql)
    data2 = c2.fetchall()
    conn2.close()

asciiout = 'synop_{m}_{p}_{ym}.txt'.format(m=m,p=p,ym=ym) 
f = open(asciiout,'w')
for l,l2 in zip(data,data2):
    if l[1]  == l2[1] and l[2]  == l2[2] and l[3]  == l2[3]:
        csv = ','.join(str(d) for d in l) + ','+str(l2[-1]); #print csv; #sys.exit(0)
        f.write(csv+'\n')
f.close()
print 'wrote ',asciiout

