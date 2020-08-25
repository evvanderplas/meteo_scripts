#!/usr/bin/env python

#  #[ imported modules
import matplotlib.pyplot as plt
import numpy
import os, sys
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.colors import colorConverter as cc
#  #]

class MySingleBox:
    def __init__(self, ax, data, xpos, label='X', width=0.5, color=None):
        #  #[ define a single box only
        self.bp = ax.boxplot(data,positions=[xpos,],widths=width)
        self.label = label
        self.xpos = xpos
        self.labelpos = xpos
        # set color
        if color:
            for key in self.bp.keys():
                clr = color
                if key=='medians': clr = 'black'
                for obj in self.bp[key]:
                    plt.setp(obj,color=clr)
        #  #]

class FilledSingleBox(MySingleBox):
    def __init__(self,*args,**kwargs):
        #  #[ define a filled box
        # extract fillcolor from kwargs
        fillcolor = kwargs.pop('fillcolor')
        # call init of parent class
        MySingleBox.__init__(self,*args,**kwargs)
        # add the fill color; code inspired on:
        # http://messymind.net/2012/07/making-matplotlib-look-like-ggplot/
        box = self.bp['boxes'][0]
        boxCoords = zip(box.get_xdata(), box.get_ydata())
        boxPolygon = Polygon(boxCoords)
        self.patch = boxPolygon
        self.patch_color = fillcolor
        #  #]

class MultiBox:
    def __init__(self, ax, data, pos=1, label='A',
                 colors=['blue','red','green','orange'], fill=False):
        #  #[ combine several boxes above one label
        num_boxes = len(data) # number of boxes per tickmark
        self.list_of_boxes = []
        
        # could be refined to use width, xspacing settings
        # and centered around input pos
        self.positions = [pos+l for l in range(num_boxes)]
        avg_pos = numpy.array(self.positions,dtype=float).mean()
        #print 'self.positions = ',self.positions
        #print 'avg_pos = ',avg_pos
        self.labelpos = avg_pos
        self.label = label
        self.colors = colors
        
        for i in range(num_boxes):
            # set alternating colors
            col_index = i % len(colors)
            if fill:
                sbp = FilledSingleBox(ax, data[i], xpos = self.positions[i],
                                      label = '', width = 0.6, 
                                      color     = colors[col_index],
                                      fillcolor = colors[col_index] )
            else:
                sbp = MySingleBox(ax, data[i], xpos = self.positions[i],
                                  label = '', width = 0.6,
                                  color = colors[col_index])

            self.list_of_boxes.append(sbp)
        #  #]
        
class MultiBoxPlot():

    def __init__(self,ax,colors = ['blue','red','green','orange']):
        self.ax = ax
        self.list_of_boxes = []
        self.labels = []
        self.label_positions = []
        self.nrofboxes = 1
        self.lastpos = -1 # if you want to have space between the boxplot groups (see also pos = self.lastpos + 2), else put to 0
        self.colors = colors
        
    def addmultibox(self, data, pos=None, label='A', fill=False):

        #if len(data) < self.nrofboxes:
        #    for em in range(self.nrofboxes - len(data)): data.append([])
        for d,dat in enumerate(data):
            #print 'data length: ',d,len(dat)#,dat
            if len(dat) == 0: 
                #print 'Empty box, moving on'
                pos = self.lastpos + 2
                self.labels.append(label)
                return 0
        
        # if no specific position is given: assume it to be next to the former, PLUS one extra to add space ...
        if pos is None: # and self.lastpos !=0: 
            pos = self.lastpos + 2
            #print 'taking former last position ',pos
        else:
            #print 'use user specified position',pos
            pass

        mb = MultiBox(self.ax, data, pos, label,
                      colors= self.colors, #['blue','red','green','orange'], 
                      fill=fill)
        self.list_of_boxes.extend(mb.list_of_boxes)
        self.labels.append(mb.label)
        self.label_positions.append(mb.labelpos)
        self.lastpos = max(mb.positions)
        if len(data) > self.nrofboxes:
            self.nrofboxes = len(data)
        if self.colors == []: # if empty take (ordered!) list
            self.colors = mb.colors
        if set(mb.colors) <= set(self.colors): # if new color, extend (ordered?) list
            self.colors.extend(set(mb.colors).difference(self.colors))

    def addbox(self, data, pos=1, label='A', fill=False, color='black'):

        # if no specific position is given: assume it to be next to the former
        if pos is None: pos = self.lastpos

        if fill:
            sbp = FilledSingleBox(self.ax, data, xpos = pos,
                                  label = label, width = 0.6, 
                                  color = color, fillcolor = color )
        else:
            sbp = MySingleBox(self.ax, data, xpos = pos,
                              label = label, width = 0.6, color = color)
        self.list_of_boxes.append(sbp)
        self.labels.append(sbp.label)
        self.label_positions.append(sbp.labelpos)

    def render(self):

        print 'where', self.label_positions

        # set ticks and labels
        self.ax.set_xticklabels(self.labels)
        self.ax.set_xticks(self.label_positions)

        #print 'Max position:', self.lastpos + 1
        self.ax.set_xlim(0,self.lastpos + 1)

        # copied from: http://matplotlib.org/examples/api/patch_collection.html
        patches = []
        patch_colors = []
        for sbp in self.list_of_boxes:
            if hasattr(sbp,'patch'):
                patches.append(sbp.patch)
                patch_colors.append(sbp.patch_color)

        if len(patches)>0:
            p = PatchCollection(patches,
                                facecolors=cc.to_rgba_array(patch_colors),
                                alpha=1.0) # 0.4)
            self.ax.add_collection(p)

        print 'where again ', self.label_positions
        return self.label_positions

    #def addlegend(self,names,loc = 1, *kwargs,**kwargs):
    def addlegend(self,names, *args,**kwargs):
        
        print 'max nr of labels: ',self.nrofboxes,self.colors
        nr    = self.nrofboxes
        nrcol = len(self.colors)
        if nrcol == 0: return 'No boxes, no legend'
        tp = {}
        for p in range(self.nrofboxes):
            #print 'wtf', p, p%nrcol, self.colors[p%nrcol], self.colors
            tp[p] = self.ax.plot([1],color = self.colors[p%nrcol],label = names[p])

        proxylegend = [tp[i][0] for i in tp.keys()]
        #self.ax.legend(proxylegend,names,loc=loc,*args,**kwargs)
        self.ax.legend(proxylegend,names,*args,**kwargs)
        

if __name__ == '__main__':
    # Some fake data to plot
    A = [[1, 2, 5,],  [7, 2], [3, 1, 2, 3, 4, 5]]
    B = [[5, 7, 2, 2, 5], [7, 2, 5, 9, 13], [2, 3, 4]]
    C = [[3, 2, 5, 7], [6, 7, 3]]
    D = [4, 2, 5, 8]
    E = [5, 4, 7, 3, 8]
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    
    pb = MultiBoxPlot(ax1)
    pb.addmultibox(A,pos=0,label='A')
    pb.addmultibox(B,pos=3,label='B',fill=True)
    pb.addmultibox(C,pos=8,label='C')
    pb.addbox(D,pos=12,label='D',color='orange')
    pb.addbox(E,pos=13,label='E',color='yellow',fill=True)
    pb.addlegend(['A','Fine','Mess'])
    pb.render()
    
    ax1.set_xlim(-2,15)
    
    filename = 'boxplot.png'
    fig.savefig(filename)
    #os.system('eog '+filename)

    print 10*'=','\ntesting loop functionality'
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    pb2 = MultiBoxPlot(ax2)
    #for pdata,lab in zip((A,B,C,D,E),(1,2,3,4,5)):
    for pdata,lab in zip((A,B,C),(1,2,3)):
        pb2.addmultibox(pdata,label=str(lab))
    print 'added boxes'
    pb2.addlegend(['diff','but','same'])
    print 'added legend'
    pb2.render()
    ax2.set_xlim(-2,15)
    filename = 'test_loop_boxplot.png'
    fig2.savefig(filename)
    #os.system('eog '+filename)

