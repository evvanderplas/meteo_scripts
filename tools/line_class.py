#!/usr/bin/env python

# 
# some
# lines 
# of 
# comment
#

#  #[ imported modules
import matplotlib.pyplot as plt
import matplotlib.cm as CM

import datetime as dt
import numpy as np
import os, sys, time


#from matplotlib.colors import colorConverter as cc
#  #]

def daterange(start_date, end_date):
    '''
    makes an interator of days within a period between startdate and enddate
    '''

    for n in range(int ((end_date - start_date).days)):
        #print start_date, dt.timedelta(days=n)
        yield start_date + dt.timedelta(days=n)


def todtg(posdate):

    #print 'posix in ',posdate
    vdate      = dt.datetime.fromtimestamp(posdate)
    isodate    = int(vdate.strftime('%Y%m%d%H')+'0000')
    
    return vdate,isodate

def toposix(dtdate):

    #print 'datetime ',dtdate,
    if type(dtdate) == dt.datetime:
        posixdate  = int(time.mktime(dtdate.timetuple()))
    
    return posixdate

class myTimeData:

    '''
    makes instance for plotting (characteristics of) different timeseries
    '''

    def __init__(self):
        self.data = {}

    def add_data(self,timedata,values):

        for t,v in zip(timedata,values):
            if not self.data.has_key(t):
                self.data[t] = []
            self.data[t].append(v)
            #if v == 0.: print 'adding',t,v,self.data[t]
                
    def dmax(self,bdate=None,edate=None):
        # gives maximum at a point in time

        times = sorted(self.data.keys())
        values = [max(self.data[t]) for t in times]
        
        return(times,values)

    def dmin(self,bdate=None,edate=None):
        times = sorted(self.data.keys())
        values = [min(self.data[t]) for t in times]
        
        return(times,values)

    def ddiff(self,bdate=None,edate=None):
        # for eg the fill_between matplotlib method
        
        times = sorted(self.data.keys())
        maxvalues = [max(self.data[t]) for t in times]
        minvalues = [min(self.data[t]) for t in times]
        
        return(times,maxvalues,minvalues)


class MySingleTimeLine:
    def __init__(self, ax, timedata, values, label='X', width=0.5, color=None):
        
        self.ax = ax
        self.label = label
        # set color
        if color is None:
            color = 'blue'
            pass # use color
        self.tsp = ax.plot(timedata,values,color = color)
        self.ax.set_xlabel('Time')


class MultiLine:
    def __init__(self, ax, timedata, values, colors=None,colpar = None,colormap = None):

        num_lines = len(values); print 'nr of lines: ', num_lines

        self.ax = ax
        if colors is not None:
            self.colors = np.array(colors)
            print 'all ',self.colors
            for l in range(num_lines):
                #print 'color ',num_lines,l,self.colors[l]
                tl = MySingleTimeLine(ax,timedata[l],values[l],color = self.colors[l])
            return None
        elif colormap is not None:
            self.colormap = CM.get_cmap(colormap)


        ## recalibrate coloring if values are not right:
        if colormap is not None and colpar is not None and len(colpar) == num_lines: # if they correspond to lines
            if np.array(colpar).dtype == '|S1': # so of type string...
                
                pass
            elif max(colpar) > 1:
                tpar = np.array(colpar,dtype=float)
                maxp = max(tpar)
                minp = min(tpar)
                print 'before',maxp,minp
                colpar = (tpar - min(tpar))/max(tpar - min(tpar))
                print 'now',max(colpar),min(colpar)
            pass
        else: # make them one color
            colpar = np.zeros(len(timedata),dtype=float) + 0.9

        #print timedata
        #print values

        for l in range(num_lines):
            #par = float(l)/num_lines; print 'color (0 < c < 256?): ',colpar[l]
            tl = MySingleTimeLine(ax,timedata[l],values[l],color = self.colormap(colpar[l]))

if __name__ == '__main__':


    postimedata = [1365415200 + 3600*lt for lt in range(0,24,3)]
    dttimedata  = [todtg(posdat)[0] for posdat in postimedata]

    postimedata2 = [1365415200 + 3600*lt for lt in range(6,36,3)]
    dttimedata2  = [todtg(posdat)[0] for posdat in postimedata2]

    valdata   = [np.sin(lt) for lt in range(0,24,3)]
    valdata2  = [np.sin(lt+1) for lt in range(6,36,3)]

    # more rigorous
    begindate = dt.datetime(2013,4,8,0)
    enddate   = dt.datetime(2013,4,10,0)
    
    # example of daterange:
    for single_date in daterange(begindate, enddate):
        print single_date.strftime('%Y-%m_%d') #dt.datetime.strftime("%Y-%m-%d", single_date)
        for st in range(0,24,3):
            sttime =  single_date + dt.timedelta(hours = st)
            print sttime

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    #lineplot = MySingleTimeLine(ax,dttimedata,valdata)
    #lineplot = MultiLine(ax,[dttimedata,dttimedata2],[valdata,valdata2],colormap = 'jet')
    lineplot = MultiLine(ax,[dttimedata,dttimedata2],[valdata,valdata2],colors = ['blue','red'])

    # pretty dates
    fig.autofmt_xdate()

    outfig = 'test_line.png'
    sf = fig.savefig(outfig)
    print 'Done',outfig
