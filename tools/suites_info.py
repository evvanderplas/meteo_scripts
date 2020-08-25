#! /usr/bin/env python

BULL_OPER_SA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'HARMONIE',
    'version':'3614',
    'model_version':'36h1.4',
    'tab':'harmonie36',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULL',
    'ltmax':48,
    }

BULL_38H11 = {'dir':'/net/bxsfsn_stc/nfs_stc/lam_oper/prodharm/scratch/hm_home/38h11/archive/{y}/{m}/{d}/{h}/',
              'shortname':'BULL_38h11',
              'fformat':'HA38_N25_{dtgh}00_{lt}00_GB',#HA38_N55_201407010600_01915_GB
              'ltmax':48,
              'tab':'harmonie37',
              'nlevs':65,
              'model':'HARMONIE',
              'version':'3811',
              'model_version':'38h1.1',
              'tab':'harmonie37',
              'name':'Harmonie (38h1.1)',
              'shortname': 'HARM38h11',
          }
BULL_OPER_LA = {
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL',
    'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'model':'HARMONIE',
    'version':'3614',
    'model_version':'36h1.4',
    'tab':'harmonie36',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULLLA',
    'ltmax':48,
    }
fourdv = {'dir':'/net/bxsfsn_ltc/nfs_ltc/lam_test/testharm/scratch/hm_home/4DVAR_300/archive/{y}/{m}/{d}/{h}/',
          'shortname':'fourdv',
          'fformat':'fc{dtgh}+{lt}grib',
          'ltmax':12,
          'tab':'harmonie37',
          'nlevs':60,
          'model':'HARMONIE',
          'version':'3712',
          'model_version':'37h1.1',
          'tab':'harmonie37',
          'name':'Harmonie (4D VAR)',
          'shortname': '4DVAR',
      }
BULL_REFO = {
    'tardir': '/data/mos/meteo/models/HARM_REFO_RESULTS',
    'tgzname': 'outLN{ymdh}.tgz', #outLN2010040200.tgz
    'dir' : '/net/bhw379/nobackup_1/users/plas/temp',
    'wdir' : '/net/bhw379/nobackup_1/users/plas/temp/{yyyy}/{mm}/{dd}/{st}',
    'fformat': 'HARM_N25_{ymdh}00_{lt}00_GB',
    'model':'HARMONIE',
    'version':'3712',
    'model_version':'37h1.2',
    'tab':'harmonie37',
    'name':'Harmonie (reforecast 37h1.2)',
    'nlevs': 60,
    'shortname': 'REFO',
    'ltmax':48,
    }

test_REFO = {
    'tardir': '/data/mos/meteo/models/HARM_REFO_RESULTS',
    'tgzname': 'outLN{ymdh}.tgz', #outLN2010040200.tgz
    'dir' : '/tmp', #'/net/bhw379/nobackup_1/users/plas/temp',
    'wdir' : '/net/bhw379/nobackup_1/users/plas/temp/{yyyy}/{mm}/{dd}/{st}',
    'fformat': 'HARM_N25_{ymdh}00_{lt}00_GB',
    'model':'HARMONIE',
    'version':'3712',
    'model_version':'37h1.2',
    'tab':'harmonie37',
    'name':'Harmonie (reforecast 37h1.2)',
    'nlevs': 60,
    'shortname': 'REFO',
    'ltmax':25,
    }

VETREFO = {
    'dir' : '/net/bhw379/nobackup_1/users/plas/temp/2010/07/14/06',
    #'wdir' : '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'HARMONIE',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETREFO',
    'ltmax':48,
    }
VETSTD = {
    'dir' : '/net/bhw379/nobackup_1/users/plas/temp/VETSTD',
    #'wdir' : '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'HARMONIE',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETSTD',
    'ltmax':18,
    }
VETP2 = {
    'dir' : '/net/bhw379/nobackup_1/users/plas/temp/VETP2',
    #'wdir' : '/net/bhw379/nobackup_1/users/plas/verif/BULL/BULL',
    #'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'HARMONIE',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'VETP2',
    'ltmax':18,
    }
