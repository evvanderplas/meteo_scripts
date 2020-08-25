#!/usr/bin/env python

import datetime as dt
today = dt.datetime.today()
td    = dt.datetime.strftime(today,'%Y%m%d')

ECJAN = {
    'dir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF',
    'wdir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF', # '/nobackup_1/users/plas/verif/BULL/ECJAN',
    'rexp': 'fc*[01][02]+*gribec',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name': 'Harmonie (ECJAN)',
    'shortname': 'ECJAN'
   }

RUC3dv = {
    'dir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/RUC3dvar',
    'wdir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/RUC3dvar', # '/nobackup_1/users/plas/verif/BULL/ECJAN',
    'rexp': 'fc*+*grib',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name': 'Harmonie (RUC)',
    'shortname': 'RUC'
   }

fourdv = {
    'dir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/4DVAR_300',
    'wdir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/4DVAR_300',
    'rexp': 'fc*+*grib',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name': 'Harmonie (4DVAR)',
    'shortname': '4DVAR',
   }

fourdv_norad = {
    'dir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/4DVAR_norad',
    'wdir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/4DVAR_norad',
    'rexp': 'fc*+*grib',
    'fformat': 'fc{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name': 'Harmonie (4DVAR)',
    'shortname': '4DVAR',
   }

MSGSIB = {
    'dir' : '/nobackup_1/users/plas/verif/BULL/MSGSIB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/MSGSIB',
    'rexp': 'ex*', #+*grib',
    'fformat': 'ex{dtgh}+{lt}grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name': 'Harmonie (MSG init)',
    'shortname': 'MSGSIB',
   }

BULL_OPER = {
    #'dir' : '/net/bhw284/nobackup/users/tijm/verg/HARM36',
    'dir' : '/data/bens03/harmonie_data/GVDB/', 
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'rexp': 'fc*'+td+'*grib_la',
    #'rexp': 'fc*[01][02]+*grib_la',
    'rexp': 'HARM_N55*[01][02]_*_GB',
    'fformat': 'HARM_N55_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL)',
    'shortname': 'BULL'
    }


BULL_OPER_SA = {
    #'dir' : '/net/bhw284/nobackup/users/tijm/verg/HARM36',
    #'dir' : '/data/bens03/harmonie_data/GVDB',
    'dir' : '/data/bens03/harmonie_data/GVDB',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'rexp': 'fc*[01][02]+*grib_sa',
    'rexp': 'HARM_N25_*[01][02]_*_GB',
    'fformat': 'HARM_N25_{dtgh}00_{lt}00_GB',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (BULL, small area)',
    'shortname': 'BULLSA'
    }

BULL_PLAS = {
    'dir' : '/nobackup_1/users/plas/verif/BULL/BULL',
    #'rexp': 'fc*'+td+'*grib_la',
    'rexp': 'ex*[01][02]+*grib',
    'model':'Harmonie',
    'model_version':'36h1.4',
    'name':'Harmonie (Oper)',
    'shortname': 'BULL'
    }

MODES = {
    'dir' : '/nobackup_1/users/plas/verif/BULL/MODES',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/MODES',
    #'rexp': 'fc*'+td+'*grib_la',
    'rexp': 'ex*[01][02]+*grib_sa',
    'model':'Harmonie',
    'version':'3614',
    'model_version':'36h1.4',
    'name':'Harmonie (Mode-S)',
    'shortname': 'MODES'
    }

BULL_37h12 = {
    'dir' : '/nobackup_1/users/plas/verif/BULL/37h12',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/37h12',
    'rexp': 'ex*grib',
    'model':'Harmonie',
    'version':'3712',
    'model_version':'37h1.2',
    'name':'Harmonie (37h1.2)',
    'shortname': '37H12'
    }

BULL_38h11 = {
    'dir' : '/net/bxsfsn_stc/nfs_stc/lam_oper/prodharm/scratch/hm_home/38h11/archive/{y}/{m}/{d}/{h}',
    'wdir' : '/nobackup_1/users/plas/verif/BULL/38h11',
    'rexp': 'HA38_N25_*00_GB',
    'fformat': 'HA38_N25_{dtgh}00_{lt}00_GB', 
    'model':'Harmonie',
    'version':'3811',
    'model_version':'38h1.1',
    'name':'Harmonie (38h1.1)',
    'shortname': '38H11'
    }


HIRLAM_D11 = {
    #'dir' : '/net/bens03/aplwsdata_cx4/GVDB',
    'dir' : '/data/bens03/aplwsdata_cx4/GVDB',
    'rexp': 'LAMH_D11*'+td+'[01][02]*', # the 00 and 12  UTC run
    'fformat': 'LAMH_D11_{dtgh}00_{lt}00_AB',
    'model':'Hirlam',
    'model_version':'7.2',
    'version':'7',
    'name':'Hirlam D11',
    'shortname': 'D11'
    }

EXP = {
    'dir' : '/net/bhw177/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF',
    'rexp': 'fc*2012030700*gribec',
    'model':'Harmonie',
    'model_version':'36h1.4',
    'name': 'Harmonie (ECJAN)',
    'shortname': 'EXP',
    'outdir' : '/nobackup/users/plas/verif/OPR/stat/plot'
    }

WIM = {
    'dir' : '/nobackup_1/users/plas/verif/BULL/WIM',
    'rexp': 'ex*',
    'fformat': 'ex{dtgh}+{lt}grib',
    'model':'Harmonie',
    'model_version':'37h1.2',
    'version':'3712',
    'name': 'Harmonie (racmo turb)',
    'shortname': 'WIM',
    }

RADAR = {
    #'dir' : '/net/bens03/aplwsdata_cx4/RADAR/forecast',
    'dir' : '/data/bens03/aplwsdata_cx4/RADAR/forecast',
    'rexp': 'RAD_NL25_PCP_FM_*00.h5', # take only hourly intensities?
    'fformat': 'RAD_NL25_PCP_FM_{dtgh}00.h5',
    'datapath':'image1/image_data',
    'metapath':'image1',
    'dtattr':'image_datetime_valid',
    'model':'Radar NL',
    'model_version':' ',
    'name': 'Dutch precipitation radar', #'European composite radar',
    'shortname': 'RADAR'
    }

RADAR_EU = {
    #'dir' : '/net/bens03/aplwsdata_cx4/RADARINPUT',
    'dir' : '/data/bens03/aplwsdata_cx4/RADARINPUT',
    'rexp': 'RAD_NL23_PCP_NA_*00.h5', # take only hourly intensities?
    'fformat': 'RAD_NL23_PCP_NA_{dtgh}00.h5',
    'datapath':'image1/image_data',
    'metapath':u'overview',
    'dtattr':u'product_datetime_end',
    'model':'Radar EUR',
    'model_version':' ',
    'name': 'European composite radar', #'European composite radar',
    'shortname': 'RADEUR'
    }

RADHDF23 = {'model':'Radar EUR',
            'name':'RAD23',
            'shortname': 'RADAR',
            #'dir':'/nobackup/users/plas/radar/', #radar_eur',
            'dir':'/nobackup_1/users/plas/verif/BULL/obs/eur', #radar_eur',
            'rexp':'RAD_NL23_PCP*h5'
            }

RADHDF25 = {'model':'Radar NL',
            'name':'RADNL_25',
            'shortname': 'RADAR',
            'dir':'/nobackup_1/users/plas/verif/BULL/obs/nl', #radar_nl',
            'rexp':'RAD_NL25_PCP*h5' # 5 min accumulated
            }

RADHDF25_3H = {'model':'Radar NL',
               'name':'RADNL_25_3H',
               'shortname': 'RADAR',
               'dir':'/nobackup/users/plas/testtemp/RAD25', 
               'rexp':'RAD_NL25_RAC*h5' # 3h accumulated
               }


MSG    = {'model':'MSG',
          'name':'MSG Cloud mask (CMa)',
          'shortname':'CM',
          'dir':'/data/mos/obsmsg/seviri/level2/safnwc/{y}/{m}/{d}/',
          'rexp':'SAFNWC_MSG2_CMa_*00_MSG*.h5',
          'fformat':'SAFNWC_MSG3_CMa__{dtgh}00_MSG-N_______.h5', #SAFNWC_MSG3_CMa__201401311000_MSG-N_______.h5
          'datapath':'CMa',
          'metapath':'',
          'dtattr':u'TIME_STAMP_UP_LINE',
          }

MSGIR  = {'model':'MSG',
          'name':'MSG IR',
          'shortname':'MSGIR',
          'dir':'/net/bhw422/nobackup/users/valkde/MSG-hdf', # /yyyy/mm/dd/
          'rexp':'METEOSAT_10_SEVIRI_EUROPE_*00_00.h5',
          'fformat':'METEOSAT_10_SEVIRI_EUROPE_{y}_{m}_{d}_{h}_00_00.h5',
          'datapath':'image5/image_data',
          'metapath':'image5/satellite',
          'dtattr':u'image_acquisition_time',
          }


VERIF   = {
           'dir' : '/nobackup/users/plas/verif/BULL/MODES',
           'rexp': 'ex*2012062112*grib',
           'model':'Harmonie',
           'model_version':'36h1.4',
           'name': 'Harmonie',
           'shortname': 'EXP',
           #'outdir' : '/nobackup/users/plas/testtemp/figs'
           'outdir' : '/nobackup/users/plas/verif/BULL/MODES/figs'
           }

ECJVER = {
    'dir' : '/net/bhw381/nobackup/users/barkmeij/harmonie/ANJAN/ECMWF',
    #'dir' : '/nobackup/users/plas/verif/BULL/ECJAN',
    'rexp': 'fc*2012062112*gribec',
    'model':'Harmonie',
    'model_version':'36h1.4',
    'name': 'Harmonie (ECJAN)',
    'shortname': 'ECJAN_VER'
   }

rad23verif = {
    'dir' : '/nobackup/users/plas/verif/BULL/obs',
    'rexp': 'RAD_NL23*00.h5', # take only hourly intensities?
    'model':'Radar EUR',
    'model_version':' ',
    'name': 'Dutch precipitation radar', #'European composite radar',
    'shortname': 'RADAR'
    }

OUT = {'dir' : './'}

sources_list = (BULL_OPER,BULL_OPER_SA,HIRLAM_D11,RUC3dv,WIM,fourdv,MSGSIB) #ECJAN to be checked!
#sources_list = (BULL_OPER,BULL_OPER_SA,HIRLAM_D11)
#sources_list = (HIRLAM_D11,)
#sources_list = (MSGSIB,WIM,)#RUC3dv,)
sources_list = (BULL_OPER,BULL_OPER_SA,HIRLAM_D11,fourdv) #,BULL_38h11)
#sources_list = (BULL_OPER,)
#sources_list = (WIM,)
#sources_list = (fourdv,RUC3dv,MSGSIB,WIM)

obs_list     = (RADAR,RADAR_EU,MSG)
obs_list     = (RADAR,RADAR_EU,MSGIR,MSG)
#obs_list     = ()
#obs_list     = (RADAR_EU,)
#obs_list     = (MSGIR,MSG,RADAR)
#obs_list     = (MSG,)
#obs_list     = (rad23verif,)


#    gribec_files = glob.glob( os.path.join(datadirjan,rexp_ECJAN))
 
