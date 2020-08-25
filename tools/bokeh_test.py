#! /usr/bin/env python

import os,sys
import numpy as np
import datetime as dt

import bokeh
import bokeh.plotting
from bokeh.objects import Range1d

import netCDF4
dtg = dt.datetime(2014,3,14,6)
dtgs = dtg.strftime("%Y%m%d%H")

ncdir      = '/net/bhw379/nobackup_1/users/plas/verif/Cabauw'
ncdir      = '/usr/people/plas/python/tools'
ncfilename = 'profile_BULLSA_{d}_De_Bilt.nc'.format(d=dtgs)
ncfile = os.path.join(ncdir,ncfilename)
print ncfile, os.path.exists(ncfile)
if not os.path.exists(ncfile): sys.exit(1)

ncobj = netCDF4.Dataset(ncfile,'r')
if 0:
    print 'In this file we have: \n',10*'='
    for k in ncobj.variables.keys(): print k, ncobj.variables[k]

time  = ncobj.variables['time'][:]
t     = ncobj.variables['T_ml'][:]
tke   = ncobj.variables['tke_ml'][:]
tcc   = ncobj.variables['tcc_ml'][:]
clice = ncobj.variables['qi_ml'][:]
clwat = ncobj.variables['ql_ml'][:]
z     = ncobj.variables['z_ml'][:]

#timeplot = np.tile([dt.datetime.fromtimestamp(t) for t in time],(60,1))
timeplot = np.tile(time,(60,1))
zplot   = z[:,::-1].T
tplot = t[:,::-1].T
tkeplot = tke[:,::-1].T
tccplot = tcc[:,::-1].T
cliplot = clice[:,::-1].T
clwplot = clwat[:,::-1].T
dts = [dt.datetime.fromtimestamp(t) for t in time]

import matplotlib.pyplot as plt
from matplotlib import cm as CM
fig = plt.figure()
ax  = fig.add_subplot(1,1,1)


cmap = CM.get_cmap('jet')
cmap = CM.get_cmap('Blues')

cllines = np.arange(1.e-4,5.1e-4,5.e-5)

#ax.contourf(timeplot,zplot,tkeplot,50)
pcont = ax.contourf(dts,z[0,::-1],tccplot,50,cmap=cmap)
pcwcont  = ax.contourf(dts,z[0,::-1],clwplot,cllines,cmap=CM.get_cmap('Reds'),alpha=0.6)
pcicont  = ax.contourf(dts,z[0,::-1],cliplot,40,cmap=CM.get_cmap('Reds'),alpha=0.6)
pclcont = ax.contour(dts,z[0,::-1],cliplot + clwplot,10,colors='red',linewidths=0.2,alpha=0.5)
tcont  = ax.contour(dts,z[0,::-1],tplot,[253,273],linewidths=2,colors='k')

ax.set_ylim(0,15000)
cb = fig.colorbar(pcwcont,shrink=0.9, extend='both',format='%.2e') #format='%.1e') 
fig.autofmt_xdate()
fig.savefig('mplprof.png')
#plt.close(fig)
sys.exit(0)

bokeh.plotting.output_file("bokprof.html", title="profile example")
bokeh.plotting.figure(tools="pan,wheel_zoom,box_zoom,reset,previewsave,select")

bokeh.plotting.image(
        image=[tkeplot], x=[0], y=[0], dw=[48], dh=[60], palette=["Spectral-11"],
        x_range = Range1d(start=0, end=48), y_range = Range1d(start=0, end=60),
        tools="pan,wheel_zoom,box_zoom,reset,previewsave", name="profile_example"
        )


if 0:
    N = 1000
    x = np.linspace(0, 4*np.pi, N)
    y = np.linspace(0, 4*np.pi, N)
    
    xx, yy = np.meshgrid(x, y)
    d = np.sin(xx)*np.cos(yy)
    
    bokeh.plotting.image(
        image=[d], x=[0], y=[0], dw=[10], dh=[10], palette=["Spectral-11"],
        x_range = Range1d(start=0, end=10), y_range = Range1d(start=0, end=10),
        tools="pan,wheel_zoom,box_zoom,reset,previewsave", name="image_example"
        )
    
#bokeh.plotting.show()  # open a browser
