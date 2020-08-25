#!/usr/bin/env python

# import the needed packages
import os, sys, glob, pickle
import time
import datetime as dt
import numpy
import matplotlib, matplotlib.pyplot
from mpl_toolkits.basemap import Basemap # map functionality
from matplotlib import cm as CM

import pygrib

def genSearchList(parameterlist):

    searchlist = []
    for parameter in parameterlist:
        if   parameter['level_indicator'] == 105: parameter['level_indicator'] = 'sfc'
        elif parameter['level_indicator'] == 109: parameter['level_indicator'] = 'ml'
        else: parameter['level_indicator'] = 'pt'

        searchlist.append([parameter['grib_indicator'],
                           parameter['level_indicator'],
                           parameter['level']])

    return parameterlist,searchlist

def translate(level_int):
    if level_int == 105:   leveltype = 'sfc'
    elif level_int == 109: leveltype = 'ml'
    elif level_int == 113: leveltype = 'pt'
    elif level_int == 8:   leveltype = 'sfc'
    elif type(level_int) == int: leveltype = 'pt' # don't know, really?
    elif type(level_int) == str: leveltype = level_int
    else: 
        print 'leveltype unknown', level_int
        leveltype = 'sfc' # default?
        #leveltype = None # marker?
        #sys.exit(1)

    return leveltype


def field_specs(param,leveltype,level,orig_values):
    
    # not used 20130131

    if param == 61 and level == 456:
        field = orig_values * 3.6e3 * 5.
    elif param == 11 and level == 2:
        field = orig_values - 273.15
    else:
        field = orig_values

    return field

def set_plotsettings(parameter,info_dict):

    #print 'set plotsettings: ',parameter

    # determine whether to use plot_levels or plot_range with colors to match
    levels = []
    if parameter.has_key('plot_levels'): # discrete list of levels
        levels = parameter['plot_levels']

    elif parameter.has_key('plot_range'):
        r = parameter['plot_range']
        try:
            levels = numpy.arange(r[0],r[1],r[2])
        except: 
            print 'Not a valid range definition: '#,r,'\nUsing default (0,10,1)'
            levels = None #numpy.arange(0,10,1)
    else:
        print 'No valid level definition: '#Using default (0,10,1)'
        levels = None #numpy.arange(0,10,1)
    

    # determine whether to use colormap or colors
    if parameter.has_key('colors'):
        colors = parameter['colors']
        cmap   = None
    elif parameter.has_key('colormap'):
        #print parameter['colormap']
        cmap = CM.get_cmap(parameter['colormap'])
        colors = None
    else: 
        colors = None
        cmap   = None 
        draw_line = True # to plot something

    # drawing lines over the contours? Default no
    if parameter.has_key('draw_line'):
        draw_line = parameter['draw_line'] 
        if parameter.has_key('line_range'):
            lr = parameter['line_range']
            try: lines = numpy.arange(lr[0],lr[1],lr[2])
            except: 
                try:
                    lines = numpy.arange(lr[0],lr[1],(lr[1]-lr[0])/len(levels))
                except:
                    lines = levels
        elif parameter.has_key('lines'):
            lines = parameter['lines']

    else: 
        draw_line = False
        lines     = None 

    # has name?
    if parameter.has_key('name'):
        name = parameter['name'] 
    else: name = 'unknown'

    # has shortname?
    if parameter.has_key('shortname'):
        shortname = parameter['shortname'] 
    else: shortname = 'unknown'

    info_dict['name']   = name
    info_dict['shortname'] = shortname
    info_dict['plotlevels'] = levels
    info_dict['colormap']   = cmap
    info_dict['plotcolors']  = colors
    info_dict['draw_line'] = draw_line
    info_dict['lines']  = lines

    return info_dict #levels,cmap,colors,draw_line,lines,name,shortname

def set_plotsettings_old(index,param_list):

    print index
    #print len(param_list), param_list[0]
    print param_list[index]

    # determine whether to use plot_levels or plot_range with colors to match
    levels = []
    if param_list[index].has_key('plot_levels'): # discrete list of levels
        levels = param_list[index]['plot_levels']

    elif param_list[index].has_key('plot_range'):
        r = param_list[index]['plot_range']
        try:
            levels = numpy.arange(r[0],r[1],r[2])
        except: 
            print 'Not a valid range definition: ',r,'\nUsing default (0,10,1)'
            levels = numpy.arange(0,10,1)
    else:
        print 'No valid level definition: Using default (0,10,1)'
        levels = numpy.arange(0,10,1)
    

    # determine whether to use colormap or colors
    if param_list[index].has_key('colors'):
        colors = param_list[index]['colors']
        cmap   = None
    elif param_list[index].has_key('colormap'):
        print param_list[index]['colormap']
        cmap = CM.get_cmap(param_list[index]['colormap'])
        colors = None
    else: 
        colors = None
        cmap   = None 
        draw_line = True # to plot something

    # drawing lines over the contours? Default no
    if param_list[index].has_key('draw_line'):
        draw_line = param_list[index]['draw_line'] 
    else: draw_line = False

    # has name?
    if param_list[index].has_key('name'):
        name = param_list[index]['name'] 
    else: name = 'unknown'

    return levels,cmap,colors,draw_line,name

def tocelsius(field):

    return field - 273.15

def post_process(grb,index,postparam_list):

    # not used? 20130131
    
    if grb.indicatorOfParameter == 33 and grb.level == 10:
        field  = numpy.sqrt( (grb.values)**2 )
        levels = postparam_list[index]['levels']
        cmap   = postparam_list[index]['cmap']
        colors = postparam_list[index]['colors']
        draw_line = postparam_list[index]['draw_line']

    elif grb.indicatorOfParameter == 61 and grb.level == 456:
        field  = (grb.values)* 3.6e3 * 5. 
        levels,cmap,colors,draw_line,name = set_plotsettings(index,postparam_list)

    elif grb.indicatorOfParameter == 11 and grb.level == 2:
        field  = (grb.values) - 273.15 
        levels,cmap,colors,draw_line,name = set_plotsettings(index,postparam_list)

    else:
        field = grb.values
        levels,cmap,colors,draw_line,name = set_plotsettings(index,postparam_list)

    return field,levels,cmap,colors,draw_line,name
