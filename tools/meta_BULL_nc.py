#! /nfs_ltc/lam_test/testharm/anaconda/bin/python

import os,sys,glob
from os.path import join as pjoin
import datetime as dt

td = dt.datetime.today()


BULL_SA  = {'dir':'/nfs_stc/lam_oper/prodharm/scratch/hm_home/36h14/archive/{y}/{m}/{d}/{st}',
         'fformat':'fc{dtgh}+{lt}grib_sa',
         'ltmax':48,
         'name':'BULL_SA',
         }
MSGSIB = {'dir':'/nfs_stc/lam_oper/prodharm/scratch/hm_home/MSG2/archive/{y}/{m}/{d}/{st}',
         'fformat':'fc{dtgh}+{lt}grib',
         'ltmax':48,
         'name':'MSGSIB',
         }
RUC   = {'dir':'/nfs_ltc/lam_test/testharm/scratch/hm_home/RUC3dvar/archive/{y}/{m}/{d}/{st}',
         'fformat':'fc{dtgh}+{lt}grib',
         'ltmax':12,
         'name':'RUC',
         }
H4DVAR = {'dir':'/nfs_ltc/lam_test/testharm/scratch/hm_home/4DVAR_300/archive/{y}/{m}/{d}/{st}',
         'fformat':'fc{dtgh}+{lt}grib',
         'ltmax':12,
         'name':'H4DVAR',
         }
WIM = {'dir':'/nfs_ltc/lam_test/testharm/scratch/hm_home/racmoturbzch02/archive/{y}/{m}/{d}/{st}',
         'fformat':'fc{dtgh}+{lt}grib',
         'ltmax':12,
         'name':'WIM',
         }

models = (BULL_SA,MSGSIB,RUC,H4DVAR,WIM)
outpath = '/nfs_ltc/lam_test/testharm/ncdat'
cmd = '/nfs_ltc/lam_test/testharm/bin/convert_harmonie.py -n -s -O -o {outfile} {inf}'

lastrun = td.replace(hour=td.hour - td.hour%3,minute=0,second=0,microsecond=0)

for ml in models:

    outdir = pjoin(outpath,ml['name'])
    try:
        os.makedirs(outdir)
    except:
        print('Outdir exists',outdir)

    for hh in range(0,24,3):
        st = lastrun - dt.timedelta(hours=hh)
        #print td,lastrun,st
        y,m,d,h = st.strftime('%Y'),st.strftime('%m'),st.strftime('%d'),st.strftime('%H')
        dtgh = st.strftime('%Y%m%d%H')
        print y,m,d,h,dtgh
        lastfile = pjoin(ml['dir'].format(y=y,m=m,d=d,st=h),ml['fformat'].format(dtgh=dtgh,lt = str(ml['ltmax']).zfill(3)))
        print lastfile, os.path.exists(lastfile)
        
        if os.path.exists(lastfile):
            files = glob.glob(pjoin(ml['dir'].format(y=y,m=m,d=d,st=h),ml['fformat'].format(dtgh=dtgh,lt = '*')))
            
            outfile = pjoin(outdir,'{name}_{dtgh}_surface.nc'.format(name=ml['name'],dtgh=dtgh))
            os.system(cmd.format(inf=files,outfile=outfile))

sys.exit(0)
bfiles = sorted(glob.glob(os.path.join(BULLPATH,'HARM_N25_*_GB')))
print bfiles[-1]
fname = os.path.split(bfiles[-1])[1]
print fname[9:19]

inf = os.path.join(BULLPATH,'HARM_N25_*'+fname[9:19]+'*_GB')
outfile = os.path.join(outpath,'BULL','harm_'+fname[9:19]+'.nc')


cmd = '/nfs_ltc/lam_test/testharm/bin/convert_harmonie.py -n -s -O -o {outfile} {inf}'.format(outfile = outfile, inf = inf)
print cmd
os.system(cmd)


