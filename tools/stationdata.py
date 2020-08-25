#!/usr/bin/env python

import os,sys,sqlite3
import numpy as np
import pandas as pd
import datetime as dt

def get_data(db,table=None,sqlin= None,begindate = None, enddate = None):

    '''
    Getting data from an SQLite database
    '''

    conn   = sqlite3.connect(db)
    c      = conn.cursor()

    info = None
    if table is not None:

        c.execute('PRAGMA table_info({tab})'.format(tab=table))
        info = c.fetchall()
    
        if sqlin == None:
            if begindate == None:
                sql = 'select * from {tab}'.format(tab = table)
            elif begindate != None:
                sql = 'select * from {tab} where '.format(tab = table)
        else:
            sql = sqlin
    elif sqlin is not None:
        sql = sqlin
    else:
        print 'What do you want from ',db
        sys.exit(1)

    try:
        c.execute(sql)
        data = c.fetchall()
    except sqlite3.OperationalError:
        data = None
    finally:
        conn.close()
    
    return info,data

def get_df(db,parameter,model,stationid = None,leadtime = None, begindate = None, enddate = None,obs=False):

    '''
    Build a query to acces sqlite HARP database, return a pandas DataFrame
    '''
    
    # resolve some inconsistencies in naming: keys is database name, values is name in obs database
    synoppar = {'apcp':'precip',
                'p':'mslp',
                't2m':'t2m',
                'wind':'wind',
                'sid':'stationid'
                }

    if type(begindate) == dt.datetime:
        import daterange
        begindate = daterange.toposix(begindate)
        enddate   = daterange.toposix(enddate)

    if obs:
        #print 'Selecting obs',db; sys.exit(0)
        sqlin   = 'select date,sdate,leadtime,stationid,{parameter} from synop'.format(parameter=synoppar[parameter],sid=str(stationid))
        columns = ['pdate','validdate','leadtime','sid',parameter+'_obs']
        stid    = 'stationid'
    else:
        sqlin = 'select date,validdate,leadtime,sid,{model} from {parameter}'.format(model=model,parameter=parameter,sid=str(stationid))
        columns= ['pdate','validdate','leadtime','sid',parameter]
        stid = 'sid'

    whr = 0 # keep track if WHERE is already used in sql

    if stationid is not None:
        sqlin += ' where {stid} = {sid}'.format(stid = stid,sid = str(stationid) )
        whr = 1 # keep track if WHERE is already used in sql

    if begindate is not None:
        if whr: sqlin += ' AND '
        else: sqlin += ' WHERE '
        sqlin += ' date BETWEEN {dtbeg} AND {dtend}'.format(dtbeg = begindate,dtend = enddate)

    if leadtime is not None:
        if whr: sqlin += 'AND'
        else: sqlin += ' WHERE '
        sqlin += ' leadtime = {lt}'.format(lt = leadtime)

        
    print 'Query database: ',db,'\nWith: ',sqlin
    info,data = get_data(db,sqlin = sqlin)
    data      = np.array(data, dtype = float) # necessary?
    print data.shape

    #dataframe = pd.DataFrame(data,columns = columns)
    try:
        return pd.DataFrame(data,columns = columns)
    except:
        print 'No data?'
        return None

def mapdatetime(dataframe):

    '''
    Takes a dataframe containing POSIX times in column 'pdate_x' (!)
    '''
    
    dataframe['datetime'] = map(dt.datetime.fromtimestamp, dataframe.pdate_x)

    return dataframe
    
def plot_scatter(dataframe,parameter,leadtime = None,
                 colorby=None,model='default',
                 outdir = './',outfile = None):

    import matplotlib.pyplot as plt
    
    # create figure canvas, axes
    fig  = plt.figure()
    ax   = fig.add_subplot(1,1,1)

    if leadtime == None:
        fc,obs = dataframe[parameter],dataframe[parameter+'_obs']
    else:
        print('Using leadtime = ',leadtime)
        fc,obs = dataframe[dataframe['leadtime_x']==leadtime][parameter],dataframe[dataframe['leadtime_x']==leadtime][parameter+'_obs']
        
    pmin,pmax = min(min(fc),min(obs)),max(max(fc),max(obs))
    if colorby is not None:
        if dataframe[colorby].shape[0] !=0:
            ax.scatter(dataframe[parameter],dataframe[parameter+'_obs'],c=dataframe[colorby])

    else:
        #ax.scatter(dataframe[parameter],dataframe[parameter+'_obs'])
        ax.scatter(obs,fc)

    
    ax.plot([pmin,pmax],[pmin,pmax],color='red')

    ax.set_title('Model {model}, parameter {par}'.format(model=model,par = parameter))
    ax.set_ylabel('forecast')
    ax.set_xlabel('observed')
    if outfile == None: 
        outfile = 'scatter_{model}_{par}.png'.format(model=model,par = parameter)
    fig.savefig(os.path.join(outdir,outfile), dpi = 150)
    print 10*'=','> created ',outfile

def plot_from_starttime(df,par = 't2m',outdir = './',outfile = None):

    import matplotlib.pyplot as plt

    # create figure canvas, axes
    fig  = plt.figure()
    ax   = fig.add_subplot(1,1,1)

    for sttime in df.datetime.unique():
        
        dfst = df[df.datetime == sttime]
        timesarr = dfst.pdate_x + dfst.leadtime_x
        times = [dt.datetime.fromtimestamp(s) for s in timesarr]
        vals  = dfst[par]
        ax.plot(times,vals)
        
    df['obstime'] = map(dt.datetime.fromtimestamp,df['pdate_y'])
    ax.scatter(df['obstime'],df[par+'_obs'],marker='+')
    ax.set_xlabel('Time')
    ax.set_ylabel(par)
    
    if outfile == None:
        outfile = 'multiseries_{p}.png'.format(p = par)
    fig.savefig(os.path.join(outdir,outfile))
    print("Created ",os.path.join(outdir,outfile))
    
    return 0

def boxplot_leadtime(df,par='t2m',outdir = './',outfile = None):

    import matplotlib.pyplot as plt
    import boxplot_class_JK as mbox

    # create figure canvas, axes
    fig  = plt.figure()
    ax   = fig.add_subplot(1,1,1)
    pb   = mbox.MultiBoxPlot(ax)

    
    data = [df[df.leadtime_x == leadtime][par] for leadtime in df.leadtime_x.unique()]
    pb.addmultibox(data,label=par)

    if outfile == None:
        outfile = 'leadtime_{p}.png'.format(p=par)
    fig.savefig(os.path.join(outdir,outfile))
    return 0


if __name__ == '__main__':

    nbk = '/nobackup/users/plas'
    datadir = '/media/plas/16GB_usb/verifdata/'
    datadir = os.path.join(os.getenv('HOME'),'HARP')
    #datadir = os.path.join(nbk,'harp')

    obsfile = os.path.join(datadir,'obs_synop_201304.db')
    fcfile  = os.path.join(datadir,'BULL_ARCH','synop_BULL_ARCH_t2m_201304.db')
    
    print 'starting up! Plotting with ipython: execfile("stationdata.py")...'
    obsthere,fcthere = os.path.exists(obsfile), os.path.exists(fcfile)
    print obsfile,obsthere,'\n',fcfile,fcthere
    if not obsthere: sys.exit(0)
    if not fcthere: sys.exit(0)

    sqlin = 'select date,sdate,leadtime,stationid,t2m from synop where stationid = 6260'
    obsinfo,obsdata = get_data(obsfile,table = 'synop',sqlin = sqlin)
    #print obsinfo
    obs = np.array(obsdata, dtype = float)
    print obs.shape

    sqlin = 'select date,validdate,leadtime,sid,BULL_ARCH from t2m where sid = 6260'
    fcinfo,fcdata = get_data(fcfile,table = 'synop',sqlin = sqlin)
    fc = np.array(fcdata,dtype = float)
    #print fcinfo

    import pandas as pd
    alldat = pd.DataFrame(fc)
    alldat.columns = ['pdate','date','leadtime','sid','t2m_fc']
    obsdat = pd.DataFrame(obs)
    obsdat[1] = obsdat[1]/10000
    obsdat.columns = ['pdate','date','leadtime','sid','t2m_ob']

    ## merge the datasets:
    #samdat = pd.merge(alldat,obsdat,how='outer') # no, 
    samdat2 = pd.merge(alldat,obsdat,on='date')
    samdat2.leadtime_x.value_counts()

    samdat2['ME'] = samdat2.t2m_fc - samdat2.t2m_ob

    figure()
    samdat2.ME[samdat2.leadtime_x == 21600].hist()
    samdat2.ME[samdat2.leadtime_x == 32400].hist(color='red')
    samdat2.ME[samdat2.leadtime_x == 75600].hist(color='green')

    samdat2['RMSE'] = np.sqrt((samdat2.t2m_fc - samdat2.t2m_ob)**2)

    lts = sorted(samdat2.leadtime_x.unique())
    somestat = []
    for lt in lts:
        mean_RMSE_lt = samdat2.RMSE[samdat2.leadtime_x == lt].mean()
        print lt,mean_RMSE_lt
        somestat.append(mean_RMSE_lt)

    import datetime as dt
    fromtimestamp = lambda x:datetime.fromtimestamp(int(x))
    samdat2['datetime'] = map(dt.datetime.fromtimestamp, samdat2.pdate_x)
    tfc  = samdat2.t2m_fc[samdat2.leadtime_x == 10800]
    tob  = samdat2.t2m_ob[samdat2.leadtime_x == 10800]
    figure()
    scatter(tob,tfc)
    plot([275,300],[275,300],color='red')

    figure()
    fcseries = pd.Series(samdat2.t2m_fc[samdat2.leadtime_x == 10800].values,samdat2['datetime'][samdat2.leadtime_x == 10800])
    obseries = pd.Series(samdat2.t2m_ob[samdat2.leadtime_x == 10800].values,samdat2['datetime'][samdat2.leadtime_x == 10800])
    fcseries.plot(color='b')
    obseries.plot(color='r')
    
