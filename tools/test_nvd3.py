#! /usr/bin/env python

import os,sys,sqlite3
from os.path import join as pjoin
import datetime as dt
import pandas as pd

import nvd3

thisdir = '/nobackup_1/users/plas/python/web/'

models = ['bull','d11','ruc','fourdv','WIM','MSGSIB']
#models = ['bull',]

outdir = '/nobackup_1/users/plas/fabriek/STAT'
topdir = '/nobackup/users/plas/harp'
modeldata   = {'bull':
              {'dir':pjoin(topdir,'BULL_3H'),
               'name':'BULL'
               },
          'd11':
              {'dir':pjoin(topdir,'D11_3H'),
               'name':'D11'
               },
          'ruc':{'dir':pjoin(topdir,'ecjan'),
                 'name':'RUC'
               },
          'fourdv':{'dir':pjoin(topdir,'ecjan'),
                 'name':'H4DVAR'
               },
          'WIM':{'dir':pjoin(topdir,'ecjan'),
                 'name':'WIM'
               },
          'MSGSIB':{'dir':pjoin(topdir,'ecjan'),
                 'name':'MSGSIB'
               },
          }
 
thresholds = [0.1,0.3,1,3]
#thresholds = [0.1,]
nbpts      = [1,3,15,49,225]
nbpts      = [15,]
           
if __name__ == '__main__':

    td = dt.datetime.today()
    yd = td - dt.timedelta(hours = 24)
    yms = td.strftime('%Y%m')
    print td.year,td.month,td.day,td.hour,yms; #sys.exit(0)

    home = os.getenv('HOME')
    db = '/home/plas/web/d3verif/cgi-bin/spatial_BULL_NL_201302.db'
    db = '/nobackup/users/plas/harp/BULL_3H/spatial_BULL_3H_201312.db'
    db = '/nobackup/users/plas/harp/BULL_3H/synop_BULL_3H_201312.db'
    
    
    p = 't2m'
    s = 6260
    m = models[0]
    
    wdir = modeldata[m]['dir']
    wnam = modeldata[m]['name']
    db = pjoin(wdir,'synop_{n}_{par}_{dtg}.db'.format(n=wnam,par=p,dtg=yms)); print db,os.path.exists(db)

    if os.path.exists(db):
        conn   = sqlite3.connect(db)
        c      = conn.cursor()
    else:
        print 'Database does not exist',db
        sys.exit(1)

    chart = nvd3.lineWithFocusChart(name='lineWithFocusChart', height=400,
                                    x_is_date=True,
                                    #x_axis_format="%_d,%m,%Y,%H%M",
                                    x_axis_format="%d %b %Y",
                                    assets_directory=thisdir)

    if p == 'fss':
        for t in thresholds:
            for nb in nbpts:
                sql = 'select * from "stats" where threshold = {thr} and nbpts = {nbp}'.format(thr=t,nbp=nb)
                c.execute(sql)
                data = c.fetchall()
        conn.close()                

    elif p in ('t2m','wind','apcp'): 
        sql = 'select * from "{p}" where sid = "{s}"'.format(p=p,s=s)
        c.execute(sql)
        data = c.fetchall()
        df = pd.DataFrame(data)
        
        print 'model',wnam
        df.columns = ['stdate','lt','validdate','sid',wnam]
        conn.close()

        for std in df.stdate.unique():
            ttimes = [dt.datetime.strptime(vd,'%Y%m%d%H') for vd in df[df.stdate == std]['validdate'].values]
            ttimes = df[df.stdate == std]['validdate'].values

            vals =  df[df.stdate == std][wnam]

            extra_serie = {"tooltip": {"y_start": "with ", "y_end": " "}, 
                               #"date_format": "%_d-%H"
                               }
            chart.add_serie(name="Parameter {p} for {m}".format(p=p,m=wnam), 
                            y=vals, x=ttimes,
                            extra=extra_serie
                            )

    outfile = os.path.join(thisdir,'d3verif','test_synop_nvd3.html')
    output_file = open(outfile, 'w')
    chart.buildhtml()
    output_file.write(chart.htmlcontent)
    output_file.close()
    print 'Created ',outfile
