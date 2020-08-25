#! /usr/bin/env python

import os,sys,sqlite3,time
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

#home = os.getenv('HOME')
#sys.path.append(os.path.join(home,'python','tools'))

def todtg(posdate):

    #print 'posix in ',posdate
    vdate      = dt.datetime.fromtimestamp(posdate)
    isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
    
    return isodate

def toposix(date):

    import time
    
    posixdate  = int(time.mktime(date.timetuple()))
    #isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
    
    return posixdate

def data_per_lt(data,lead,ind_lt,ind_par):

    dlt = [d[ind_par] for d in data if d[ind_lt] == lead]

    return dlt

def data_per_nbp(data,nbp,ind_nbp,ind_par):

    dnbp = [d[ind_par] for d in data if d[ind_nbp] == nbp]

    return dnbp

def tovaliddate(date,ltmax=24,lts=3):
    '''
    creates a list of startingdates and leadtimes that are valid at the same time
    eg if date is dt.datetime(2013,7,2,6) it gives ['2013070106+024','2013070109+021','2013070112+018',..]
    to use in SQL "where validdate in (...)" or so
    '''
    
    vdate = []
    for l in range(lts,ltmax,lts):
        stdate = date - dt.timedelta(hours = l)
        dtg,lt = stdate.strftime('%Y%m%d%H'),str(l).zfill(3)
        vdate.append('"'+dtg+'+'+lt+'"')
    
    return vdate


def get_data(db,table,sqlin= None,begindate = None, enddate = None,verb=False):

    '''
    Getting data from an SQLite database
    '''

    conn   = sqlite3.connect(db)
    c      = conn.cursor()

    c.execute('PRAGMA table_info({tab})'.format(tab=table))
    info = c.fetchall()
    
    if sqlin == None:
        if begindate == None:
            sql = 'select * from {tab}'.format(tab = table)
        elif begindate != None:
            sql = 'select * from {tab} where '.format(tab = table)
    else:
        sql = sqlin
            
    if verb: print 'From ',db,'execute ',sql
    try:
        c.execute(sql)
        data = c.fetchall()
    except:
        print 'Error '
        return None,None
    finally:
        conn.close()
    
    return info,data

def plotdict(datadict,labels = None,outfig = './somefig.png',title = 'A boxplot from a dictionary'):

    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    
    ax.set_title(title)

    #sys.path.append('/usr/people/plas/python/tools')
    import boxplot_class_JK as mbox

    if labels == None: labels = sorted(datadict.keys())
    if len(labels) < len(datadict.keys()):
        print 'Warning: too few labels, using keys, \nlabels:',labels,'\nkeys: ',datadict.keys()
        labels = sorted(datadict.keys())

    pb = mbox.MultiBoxPlot(ax)
    for f,lab in zip(sorted(datadict.keys()),labels):
         pb.addmultibox([datadict[f]],label=str(lab))
    pb.render()

    plotfile = outfig
    sf = fig.savefig(plotfile)
    print 'plotted ',plotfile
    plt.close(fig)
    
def plotmultidict(datadict,labels = None,outfig = './somefig.png',title = 'Grouped box plot'):

    '''
    makes a boxplot of different datasets, expects a dict 
    '''


    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    ax.set_title(title)

    #sys.path.append('/usr/people/plas/python/tools')
    import boxplot_class_JK as mbox

    if labels == None: labels = datadict.keys()
    if len(labels) < len(datadict.keys()):
        print 'Warning: too few labels, using keys, '+'\nlabels:',labels,'\nkeys: ',datadict.keys()
        labels = datadict.keys()

    pb = mbox.MultiBoxPlot(ax)

    cases = sorted(datadict.keys()) ## NB cases are sorted alphabetically!
    if str in [type(c) for c in cases]:
        print 'multidict: ',cases,datadict[cases[0]].keys()
        labels = sorted(datadict[cases[0]].keys())

        # take the first case as determining the labels etc(?)
        for f,lab in zip(sorted(datadict[cases[0]].keys()),labels):
            # [datadict[c1][f],datadict[c2][f],...]
            pb.addmultibox([datadict[c][f] for c in cases],label=str(lab))

        pb.addlegend(cases,loc=4)
            
    else: # assume it is a 1 case dict
        print 'one case dict: ',datadict.keys()
        for f,lab in zip(sorted(datadict[c].keys()),labels):
            pb.addmultibox([datadict[f]],label=str(lab))

    pb.render()

    plotfile = outfig
    sf = fig.savefig(plotfile)
    print 'plotted ',plotfile
    plt.close(fig)
    
def boxplotNoDict(dbfiles,cases,thresholds,nbpts,param = 'fss',tab = 'stats',colors=None,baserateThreshold = 0.2,outdir = './'):

    '''
    makes comparative boxplots of param (default 'fss') based on the
    HARP spatial databases "dbfiles" (eg 'spatial_BULL_3H_201310.db'),
    described by "cases" (eg 'BULL','D11' etc), where the results are
    generated with given thresholds and neighbourhoods (nbpts).

    The table where this is set is default set to 'stats'.

    A baserateThreshold (default set to 0.2, so 20% of the pixels has
    rain > threshold) is used to filter out so-called 'non-events' as
    this gives a lot of points with value 0 in the distribution.
    '''

    try:
        os.makedirs(outdir)
    except:
        print 'output to ',outdir

    if colors == None: colors = ['blue','red','green','purple','orange']

    sys.path.append('/usr/people/plas/python/tools')
    import boxplot_class_JK as mbox

    print cases
    print dbfiles
    #sys.exit(1)
 
    outfiles = []
    for thr in thresholds:

        fignb = plt.figure()
        axnb  = fignb.add_subplot(1,1,1)
        pbnb  = mbox.MultiBoxPlot(axnb, colors = colors)
        axnb.set_title('{PAR} (> {THR} mm/3h) as a function of NBR SIZE'.format(PAR = param, THR = thr))

        for nbp in nbpts:

            ## nbp plot: eg FSS as a function of neighbourhood size
            datanbp = []
            for c,db in zip(cases,dbfiles):
                sql_fss_nbp = 'select {PAR} from {TAB} '+'where threshold = {THR} and nbpts = {NBP} and baserate > {BR}'
                sqlin = sql_fss_nbp.format(PAR = param, TAB = tab, 
                                           THR=str(thr), NBP = str(nbp), 
                                           BR = str(baserateThreshold))
                print sqlin
                inf,sdat = get_data(db,tab,sqlin = sql_fss_nbp.format(PAR = param, TAB = tab, 
                                                                      THR=str(thr), NBP = str(nbp), 
                                                                      BR = str(baserateThreshold)))
                datanbp.append(np.array([d[0] for d in sdat])) # append to make new list within a list

            pbnb.addmultibox(datanbp,label=str(nbp))

            ## leadtime plot: eg FSS as a function of leadtime
            figlt = plt.figure()
            axlt  = figlt.add_subplot(1,1,1)
            pblt  = mbox.MultiBoxPlot(axlt, colors = colors)
            axlt.set_title('{PAR} (> {THR} mm/3h) as a function of leadtime, size {NBPTS}'.format(PAR = param, THR = thr, NBPTS=nbp))

            for lt in range(3,49,3):

                datalt = []
                for c,db in zip(cases,dbfiles):

                    sql_gen =  'select fss from stats '+\
                        'where threshold = {THR} and nbpts = {NBP} '+\
                        'and leadtime = {LT} and baserate > {BR}'

                    lts = 3600*lt # POSIX: from hours to seconds

                    # getting the actual data
                    sql_fss_lt = sql_gen.format(THR=str(thr), NBP = str(nbp), LT=str(lts), BR = str(baserateThreshold))                   
                    inf,sdat = get_data(db,tab,sqlin = sql_fss_lt)
                    datalt.append(np.array([d[0] for d in sdat]))
                    print 'Case:',c,db,sql_fss_lt

                print 'for leadtime {LT} data: {DAT}'.format(LT = str(lt),DAT = datalt)
                pblt.addmultibox(datalt,label=str(lt))

            pblt.addlegend(cases,loc=3)
            pblt.render()

            outfig = os.path.join(outdir,'{PAR}_thr{THR}_nbp{NBP}_lt.png'.format(PAR=param,THR=thr,NBP=nbp))
            plotfile = outfig
            sf = figlt.savefig(plotfile)
            print '+++++++++ plotted ',plotfile
            plt.close(figlt)
            outfiles.append(plotfile)

        # end loop over NBPTS
        pbnb.addlegend(cases,loc=4)
        pbnb.render()

        outfig = os.path.join(outdir,'{PAR}_thr{THR}_nb.png'.format(PAR=param,THR=thr))
        plotfile = outfig
        sf = fignb.savefig(plotfile)
        print '+++++++++ plotted ',plotfile
        outfiles.append(plotfile)
        plt.close(fignb)
        
    # end loop over thresholds
    return outfiles
    
def matrixplot(dbfiles,cases,thresholds,nbpts,param = 'fss',tab = 'stats',baserateThreshold = 0.2,
                     begindate = None, enddate = None,
                     colormap = None,
                     outdir = './'):

    if len(dbfiles) != len(cases): 
        print 'One case per database, presumably?'
        print dbfiles, cases
        sys.exit(1)

    try:
        os.makedirs(outdir)
    except:
        print 'output to ',outdir


    import matplotlib.cm as cm   

    sql_gen =  'select {PAR} from {TAB} '+\
        'where threshold = {THR} and nbpts = {NBP} '+\
        'and baserate > {BR}'
    if begindate is not None and type(begindate) == dt.datetime:
        sql_gen = sql_gen + ' and date > {bd}'.format(bd = toposix(begindate))
    if enddate is not None and type(enddate) == dt.datetime:
        sql_gen = sql_gen + ' and date < {ed}'.format(ed = toposix(enddate))
    
    outfiles = []
    for db,c in zip(dbfiles,cases):
        
        X = np.zeros((len(nbpts),len(thresholds)),dtype=float)
        fig,ax = plt.subplots()
        for i,thr in enumerate(thresholds):
            
            for j,nbp in enumerate(nbpts):
            
                sql = sql_gen.format(PAR=param,TAB=tab,THR=thr,NBP=nbp,BR=baserateThreshold)
                print 'Query: ',db,sql
                inf,sdat = get_data(db,tab,sql)
                print c,thr,nbp,len(sdat), type(sdat)
                X[j,i] = np.array(sdat).mean()
                print X

        ax.imshow(X, cmap=cm.jet, interpolation='nearest',vmin=0.25,vmax=1.,alpha=0.7)
        for i in range(X.shape[1]): # thresholds
            for j in range(X.shape[0]): # nbpts
                plt.text(i,j, '{0:1.2f}'.format(X[j,i]),ha='center',va='center')

        ax.set_xticks(range(X.shape[1]))
        ax.set_xticklabels(thresholds)
        ax.set_xlabel('Threshold')

        ax.set_yticks(range(X.shape[0]))
        ax.set_yticklabels(nbpts)
        ax.set_ylabel('Neighbourhood')
        ax.set_title('{PAR} for {case} as a function of \nthreshold, neighbourhood'.format(PAR=param,case=c))
        outfile = os.path.join(outdir,'matrix_{case}_{ymd}.png'.format(case=c,ymd = begindate.strftime('%Y%m%d')))
        fig.savefig(outfile)
        print 10*'=',"> Created ",outfile
        plt.close(fig)
        outfiles.append(outfile)
    return outfiles


def lineplot_spatial(dbfiles,cases,thresholds,nbpts,param = 'fss',tab = 'stats',baserateThreshold = 0.2,
                     begindate = None, enddate = None,
                     colors = None,
                     colormap = None,
                     outdir = './', verb= False):

    '''
    Makes an array of timeseries plots for eg fraction skill scores a
    s a function of time.  Warning: to make it easy on the eye, for
    each point the maximum is taken, and a fill plot is used.
    
    '''

    ## plots that do not appeal very much
    per_case = False
    single_t = False

    try:
        os.makedirs(outdir)
    except:
        print 'output to ',outdir

    if colors == None: 
        colors = ['blue','red','darkgreen','purple','green','orange'] 
        #colors = ['blue','red','green','purple','orange']


    import pandas as pd


    # init plots
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    ax2 = ax.twinx() # for baserate bars!
    ax.set_ylim(0,1.1)
    ax2.set_ylim(0,1.1)

    #sys.path.append('/usr/people/plas/python/tools')
    #import line_class as mline

    print cases
    print dbfiles

    sql_gen =  'select date,leadtime,{PAR},baserate from {TAB} '+\
        'where threshold = {THR} and nbpts = {NBP} '+\
        'and date = {STDATE} and baserate > {BR}'
    
    sql_gen =  'select date,leadtime,{PAR},baserate from {TAB} '+\
        'where threshold = {THR} and nbpts = {NBP} '+\
        'and validdate in ({VDATE}) and baserate > {BR}'

    lts   = 3
    ltmax = 24

    if begindate is not None:
        
        dtg = begindate.strftime('%Y%m%d')

        if enddate is not None:
            rng = pd.date_range(begindate,enddate,freq='3H')
            plotdict = {}

            # keep track of produced plots
            outfiles = []

            # loop over nbpts! 
            for nbp in nbpts:

                figmt  = plt.figure()
                axmt   = {}

                for t,thr in enumerate(thresholds):

                    if single_t:
                        # plot for one threshold, compare models
                        figm = plt.figure()
                        axm   = figm.add_subplot(2,1,1)
                        axmf  = figm.add_subplot(2,1,2)

                    # plot for one threshold, compare models
                    axmt[t] = figmt.add_subplot(len(thresholds),1,t+1)

                    plotdict[thr] = {}
                    for c,db,col in zip(cases,dbfiles,colors):
                    
                        tdat,brdat = [],[] # temporary list
                        for r in rng: 
                            # print r
                    
                            vdate = tovaliddate(r,ltmax = ltmax,lts = lts)
                            
                            ## .format(**mydict)
                            sqlline = sql_gen.format(PAR=param,TAB=tab,THR=thr,NBP=nbp,BR=baserateThreshold,VDATE=','.join(vdate))
                            #print sqlline

                            # get date,leadtime, parameter and baserate
                            inf,sdat = get_data(db,tab,sqlin = sqlline)

                            if sdat != []:
                                fssdat   = [s[2] for s in sdat]
                                baserate = [s[3] for s in sdat]
                                if verb: print r, fssdat, max(fssdat), baserate
                                tdat.append(max(fssdat))
                                brdat.append(max(baserate))
                            else:
                                tdat.append(0.)
                                brdat.append(0.)

                        plotdict[thr][c] = pd.Series(tdat,index=rng)
                        #print c,plotdict[thr][c]

                        if single_t:
                            # plot the line in model comparison plot
                            plotdict[thr][c].plot(ax=axm,color=col,label=c)

                            axmf.fill_between(rng,tdat,alpha=0.4,color=col)
                            axmf.fill_between(rng,brdat,alpha=0.1,color='gray')

                        print 'plotting ',t,' for case ',c,col
                        axmt[t].set_ylim(0,1)
                        axmt[t].fill_between(rng,tdat,alpha=0.4,color=col,label=c)
                        axmt[t].plot(rng,tdat,alpha=0.4,color=col,label=c)
                        axmt[t].fill_between(rng,brdat,alpha=0.1,color='gray')
                        axmt[t].set_ylabel('pcp > {thr}'.format(thr=thr))
                    

                    if single_t:
                        axm.legend(loc='best')
                        outfile = os.path.join(outdir,'lineplot_cmpmodels_{t}_{nb}_{d}.png'.format(t=thr,nb=nbp,d=dtg))
                        figm.savefig(outfile)
                        print 10*'=','> created ',outfile
                        plt.close(figm)
                        
                #plt.setp([a.get_xticklabels() for a in axmt[:-1]], visible=False)
                #plt.setp(axmt[0].get_xticklabels(), visible=False)

                axmt[3].legend(loc='best',fontsize=9)
                figmt.suptitle('FSS for thresholds {thrs}, nb size {nb} (at {dtg}) '.format(thrs = thresholds,nb=nbp,dtg=rng[-1].strftime('%Y%m%d')))
                figmt.autofmt_xdate()
                outfile = os.path.join(outdir,'lineplot_cmpmodels_mul_nb{nb}_{d}.png'.format(nb=nbp,d=dtg))
                figmt.savefig(outfile)
                print 10*'=','> created ',outfile; #sys.exit(1)
                plt.close(figmt)
                outfiles.append(outfile)

                
                if per_case:
                    # plot for one model, compare thresholds
                    for c in cases:
                        figt = plt.figure()
                        axt  = figt.add_subplot(1,1,1)
                        
                        for thr,col in zip(thresholds,colors):
                            # plot the line in threshold comparison plot
                            plotdict[thr][c].plot(ax=axt,label=str(thr),color=col)
                        axt.legend(loc='best')

                        outfile = os.path.join(outdir,'lineplot_cmpthresholds_{c}_{nb}_{d}.png'.format(c=c,nb=nbp,d=dtg))
                        figt.savefig(outfile)
                        print 10*'=','> created ',outfile
                        plt.close(figt)
                        
        #sys.exit(1)
    return outfiles
    


def lineplotNoDict(dbfiles,cases,thresholds,nbpts,param = 'fss',tab = 'stats',
                   baserateThreshold = 0.2,
                   begindate = None, enddate = None,
                   colormap = None,
                   outdir = './'):

    try:
        os.makedirs(outdir)
    except:
        print 'output to ',outdir


    import matplotlib.pyplot as plt
    import numpy as np
    sys.path.append('/usr/people/plas/python/tools')
    import line_class as mline

    print cases
    print dbfiles
    #sys.exit(1)

    for thr in thresholds:

        for nbp in nbpts:

            # init plots
            fig = plt.figure()
            ax  = fig.add_subplot(1,1,1)
            ax2 = ax.twinx() # for baserate bars!
            ax.set_ylim(0,1.1)
            ax2.set_ylim(0,1.1)


            figt = plt.figure()
            axt  = figt.add_subplot(1,1,1)
            axt.set_ylabel(param)
            axt.set_ylim(-0.1,1.1)

            figfb = plt.figure()
            axfb  = figfb.add_subplot(1,1,1)
            axfb.set_ylabel(param)
            axfb.set_ylim(-0.1,1.1)

            tt1 = mline.myTimeData()
            tt2 = mline.myTimeData()
            tt3 = mline.myTimeData()

            timedat,valdat,col = [],[],[]
            for stdate in mline.daterange(begindate, enddate):

                for sthour in range(0,24,3):

                    stdt = stdate + dt.timedelta(hours = sthour)
                    st   = mline.toposix(stdt)

                    sql_gen =  'select date,leadtime,{PAR},baserate from {TAB} '+\
                        'where threshold = {THR} and nbpts = {NBP} '+\
                        'and date = {STDATE} and baserate > {BR}'

                    sqlline = sql_gen.format(PAR=param,TAB=tab,THR=thr,NBP=nbp,BR=baserateThreshold,STDATE=st)

                    colors = ['blue','red','green'] #['green','blue','red']
                    colors = ['blue','red','darkgreen','purple','green','orange']
                    #for c,db,ccolor in zip(cases,dbfiles,colors): 
                    for c,db,ccolor,tt in zip(cases,dbfiles,colors,(tt1,tt2,tt3)): 

                        inf,sdat = get_data(db,tab,sqlin = sqlline)

                        times  = np.array([mline.todtg(d[0]+d[1])[0] for d in sorted(sdat)])
                        values = np.array([d[2] for d in sorted(sdat)])
                        baser  = np.array([d[3] for d in sorted(sdat)])

                        # col.extend([st])  # if coloring by starting time
                        # col.extend([nbp]) # if coloring by neighbourhood
                        # col.extend([thr]) # if coloring by threshold
                        col.extend([ccolor])                       
                        
                        # good idea? forceFill fills in a non-existent value with 0.
                        forceFill = True
                        if len(values) == 0 and forceFill:
                            #print stdt, 0.
                            tt.add_data([stdt],[0.])
                        else:
                            tt.add_data(times,values)
                        
                        ax.plot(times,values,color = ccolor,label=c)

                    if times != []:
                        #print times,baser
                        ax2.fill_between(times,baser,0,color='gray',alpha=0.3)

            ## end loop over starting dates
            
            for tt,col,c in zip((tt1,tt2,tt3),colors,cases):
                axt.plot(*tt.dmax(),color = col,alpha = 0.9,label = c) 
                axt.fill_between(*tt.dmax(),color = col,alpha = 0.4) 
                axfb.plot(*tt.dmax(),color = col,alpha = 0.9,label = c) 
                axfb.fill_between(*tt.ddiff(),color = col,alpha = 0.4) 

            figt.autofmt_xdate()
            axt.legend(loc='best')
            axt.set_title('{PAR} as a function of time, thr = {THR}, nbp = {NBP}'.format(PAR=param, THR=thr,NBP=nbp))

            figfb.autofmt_xdate()
            axfb.legend()
            axfb.set_title('{PAR} as a function of time, thr = {THR}, nbp = {NBP}'.format(PAR=param, THR=thr,NBP=nbp))

            #axt.fill_between(*[tt1.dmax(),tt2.dmax(),tt3.dmax()],color = colors,alpha = 0.4) 
            fbfile = os.path.join(outdir,'time_{PAR}_thr{THR}_nbp{NBP}_timefb.png'.format(PAR=param,THR=thr,NBP=nbp))
            figfb.savefig(fbfile)
            plt.close(figfb)
            maxfile = os.path.join(outdir,'time_{PAR}_thr{THR}_nbp{NBP}_timemax.png'.format(PAR=param,THR=thr,NBP=nbp))
            figt.savefig(maxfile)
            print 'plotted testdata',thr,nbp
            plt.close(figt)
            #sys.exit(0)

            # format the ticks
            import matplotlib
            if 0:

                # pretty
                years    = matplotlib.dates.YearLocator()   # every year
                months   = matplotlib.dates.MonthLocator()  # every month
                days     = matplotlib.dates.DayLocator()    # every day
                lhours   = matplotlib.dates.HourLocator(np.arange(6,25,6))   # every hour
                yearsFmt = matplotlib.dates.DateFormatter('%Y')
                daysFmt  = matplotlib.dates.DateFormatter('%Y-%m-%d') #, %H:%M')
                hoursFmt = matplotlib.dates.DateFormatter('%H')
                
                matplotlib.pyplot.setp(ax.get_xticklabels(), visible=False)
                
                ax.xaxis.set_major_locator(days)
                ax.xaxis.set_major_formatter(daysFmt)
                ax.xaxis.set_minor_locator(lhours)
                ax.xaxis.set_minor_formatter(hoursFmt)

            #days     = matplotlib.dates.DayLocator()    # every day
            #daysFmt  = matplotlib.dates.DateFormatter('%Y-%m-%d') #, %H:%M')

            fig.autofmt_xdate()
            ax.legend()

            outfig = os.path.join(outdir,'{PAR}_thr{THR}_nbp{NBP}_time.png'.format(PAR=param,THR=thr,NBP=nbp))
            plotfile = outfig
            sf = fig.savefig(plotfile)
            print '+++++++++ plotted ',plotfile
            plt.close(fig)
            #sys.exit(1)


for i in range(1):
    print i

    
if __name__ == '__main__':

    if os.path.isdir('/nobackup/users/plas'):
        harpdir = '/nobackup/users/plas/harp'
    else:
        harpdir = '/home/plas/HARP'

    fcdb = os.path.join(harpdir,'BULL_ARCH/spatial_BULL_3H_201305.db')

    fcdb_BULL  = os.path.join(harpdir,'BULL_ARCH','spatial_BULL_3H_201304.db')
    fcdb_MODES = os.path.join(harpdir,'MODES_3H','spatial_MODES_201304.db')
    fcdb_D11   = os.path.join(harpdir,'D11_ARCH','spatial_D11_3H_201304.db')

    tab  = 'stats'
    
    print 'Getting data ',fcdb,tab
    inf,data = get_data(fcdb,tab)
    for col in inf: print col[0],col[1]
    names = [ n[1] for n in inf]
    dtg = names.index('date')
    lt  = names.index('leadtime')
    fss = names.index('fss')
    baser = names.index('baserate')

    print dtg,lt,data[0][dtg],data[0][lt]
    print type(data),len(data),len(data[0])

    for d in data[0:10]:
        print d[dtg],d[lt],todtg(d[dtg]+d[lt])

    posdtgex = 1366797600

    


    print 'selecting specific records, thr=0.1, nbp=3'
    sqlprefab = 'select * from stats where threshold = 0.1 and nbpts = 3 and baserate > 0.1'
    inf_B,data_B = get_data(fcdb_BULL,tab,sqlin = sqlprefab)
    inf_M,data_M = get_data(fcdb_MODES,tab,sqlin = sqlprefab)
    inf_D,data_D = get_data(fcdb_D11,tab,sqlin = sqlprefab)

    for d in data_B[0:10]:
        print d[dtg],d[lt],todtg(d[dtg]+d[lt]),d[fss]

    fssdat   = [data_per_lt(data_B,3600*lead,lt,fss) for lead in range(3,49,3)]
    baserdat = [data_per_lt(data_B,3600*lead,lt,baser) for lead in range(3,49,3)]

    fssdat_M   = [data_per_lt(data_M,3600*lead,lt,fss) for lead in range(3,49,3)]
    baserdat_M = [data_per_lt(data_M,3600*lead,lt,baser) for lead in range(3,49,3)]

    fssdat_D   = [data_per_lt(data_D,3600*lead,lt,fss) for lead in range(3,49,3)]

    #print fssdat
    #sys.exit(0)
    
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax  = fig.add_subplot(2,1,1)

    bp = matplotlib.pyplot.boxplot(fssdat, notch=0, sym='+', vert=1, whis=1.5)

    ax  = fig.add_subplot(2,1,2)

    #home = os.getenv('HOME')
    #sys.path.append(os.path.join(home,'python','tools')
    import boxplot_class_JK as mbox

    pb = mbox.MultiBoxPlot(ax)
    for f,b,d,lt in zip(fssdat,fssdat_M,fssdat_D,range(3,49,3)):
        print len(f),len(b)
        if len(f) == 0 or len(b) == 0 or len(d) == 0: break
        pb.addmultibox([f,b,d],label=str(lt))
        
    #pb.addlegend(['FSS','Baserate'])
    pb.addlegend(['FSS BULL','FSS Mode-S','FSS D11'])
    pb.render()

    #bp = matplotlib.pyplot.boxplot(fssdat+baserdat, notch=0, sym='+', vert=1, whis=1.5)
    
    plotfile = 'box_pyclass.png'
    sf = fig.savefig(plotfile)
