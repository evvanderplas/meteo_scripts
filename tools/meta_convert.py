#! /usr/bin/env  python

import os,sys,glob
import datetime as dt
import subprocess
import gc

sys.path.append('/usr/people/plas/python/tools')
import nc_class as ncprof

td = dt.datetime.today()

#BULLPATH = '/net/bens03/harmonie_data/GVDB/'
suites = {
    'BULL_OPER_SA':{'dir':'/data/bens03/harmonie_data/GVDB/',
            'fformat':'HARM_N25_{dtgh}00_{lt}00_GB',
            'ltmax':48,
            'tab':'harmonie36',
            'nlevs':60,
            'model':'Harmonie',
            'version':'3614',
            'model_version':'36h1.4',
            'name':'Harmonie (BULL, small area)',
            'shortname': 'BULL',
        },
    'BULL_38H11':{'dir':'/net/bxsfsn_stc/nfs_stc/lam_oper/prodharm/scratch/hm_home/38h11/archive/{y}/{m}/{d}/{h}/',
                  'shortname':'BULL_38h11',
                  'fformat':'HA38_N25_{dtgh}00_{lt}00_GB',#HA38_N55_201407010600_01915_GB
                  'ltmax':48,
                  'tab':'harmonie37',
                  'nlevs':65,
                  'model':'Harmonie',
                  'version':'3811',
                  'model_version':'38h1.1',
                  'name':'Harmonie (38h1.1)',
                  #'shortname': 'HARM38h11',
              },
    'fourdv':{
        #'dir':'/net/bxsfsn03/nfs_ltc/lam_test/testharm/scratch/hm_home/4DVAR_300/archive/{y}/{m}/{d}/{h}/',
        'dir':'/net/bxsfsn_ltc/nfs_ltc/lam_test/testharm/scratch/hm_home/4DVAR_300/archive/{y}/{m}/{d}/{h}/',
        'shortname':'fourdv',
        'fformat':'fc{dtgh}+{lt}grib',
        'ltmax':12,
        'tab':'harmonie37',
        'nlevs':60,
        'model':'Harmonie',
        'version':'3712',
        'model_version':'37h1.1',
        'name':'Harmonie (4D VAR)',
        'shortname': '4DVAR',
        },
    }

#ECJANPATH = '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF'
models = ['BULL_OPER_SA','BULL_38H11','fourdv']
#models = ['BULL_38H11',]
#models = ['fourdv',]

outpath = '/net/bhw379/nobackup_1/users/plas/ncdat'
outpath = '/nobackup_1/users/plas/ncdat'


today = dt.datetime.today()
today = today.replace(hour=today.hour - today.hour%3,minute=0,second=0,microsecond=0)
td = dt.datetime.strftime(today,'%Y%m%d%H')

from memory_inspector import report_mem_usage
memtest = True
if memtest: report_mem_usage() # init

for mod in models:

    print mod
    for st in range(0,30,3):

        starttime = today - dt.timedelta(hours=st)

        y,m,d,h   = starttime.strftime('%Y'),starttime.strftime('%m'),starttime.strftime('%d'),starttime.strftime('%H')
        sdtg      = starttime.strftime('%Y%m%d%H')

        sdir = suites[mod]['dir'].format(y=y,m=m,d=d,h=h)
        lastfile = suites[mod]['fformat'].format(dtgh = sdtg,lt = str(suites[mod]['ltmax']).zfill(3))

        print 'Exists:',os.path.join(sdir,lastfile),os.path.exists(os.path.join(sdir,lastfile))
        if not os.path.exists(os.path.join(sdir,lastfile)): continue
    
        
        sfc = False
        sfc = True
        if sfc:

            # check if already done:
            outfile = os.path.join(outpath,suites[mod]['shortname'],'harm_{model}_{dtgh}.nc'.format(model=suites[mod]['shortname'],dtgh = sdtg))
            print outfile, os.path.exists(outfile)
            if not os.path.exists(outfile): 

                if memtest:
                    print 'Surface loop over leadtimes'
                    report_mem_usage() # report

                inf = os.path.join(sdir,suites[mod]['fformat'].format(dtgh = sdtg,lt='*'))

                cmd = '/usr/people/plas/python/tools/convert_harmonie_ext.py -n -s -l {levs} -t {tab} -prt -O -o {outfile} {inf}'.format(tab=suites[mod]['tab'], outfile = outfile, inf = inf,levs=suites[mod]['nlevs'])
                print cmd

                # run command
                pid = subprocess.Popen(cmd,shell=True)
                pid.communicate()
                #os.system(cmd)

        if memtest:
            print 'Done: Surface loop over leadtimes'
            report_mem_usage() # report

        
        prof = True
        #prof = False
        if prof:

            if memtest:
                print 'Profile loop over leadtimes'
                report_mem_usage() # report


            cmd = 'python /usr/people/plas/python/tools/use_ncprof_cli_station.py {model} {startdate} {indir} {outdir}'.format(model = mod, startdate=starttime.strftime('%Y%m%d%H'),indir = sdir, outdir = os.path.join(outpath,suites[mod]['shortname'])) 
            pid = subprocess.Popen(cmd,shell=True)
            pid.communicate()

            if memtest:
                print 'Done: Profile loop over leadtimes'
                report_mem_usage() # report


sys.exit(0)

