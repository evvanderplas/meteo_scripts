def do_kdtree(combined_x_y_arrays,points):
    import scipy.spatial
    
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes

def toLatlon(latlons,index):

    return [(latlons[0][i],latlons[1][i]) for i in index]
    #return [(latlons[0][i],latlons[1][i]) for (i[1],i[0]) in index]

def find_four_neighbours(latlons,(xlat,xlon)):

    lats,lons = latlons[0],latlons[1]
    flats,flons = latlons[0].ravel(),latlons[1].ravel()
    dim0,dim1   = latlons[0].shape

    llarray = numpy.dstack([flats,flons])[0]
    res2    = do_kdtree(llarray,(xlat,xlon))
    i0,j0 = res2%dim0, (res2-res2%dim0)/dim0
    print 'KDTree: ',i0,j0,flats[res2],flons[res2]

    from matplotlib.path import Path
    
    jn0,in0 = res2%dim0, (res2-res2%dim0)/dim0
    sq1   = [(in0,jn0),( in0+1,jn0),(in0+1,jn0+1),(in0,jn0+1)]
    print sq1
    llsq1 = toLatlon(latlons,sq1)
    print llsq1
    path1 = Path(llsq1)
    print 'First test: ',path1,path1.contains_point((xlat,xlon))
    sq = {}
    sq['ll']   = [(in0,jn0),(in0+1,jn0),(in0+1,jn0+1),(in0,jn0+1)]
    sq['lr']   = [(in0-1,jn0),(in0,jn0),(in0,jn0+1),(in0-1,jn0+1)]
    sq['ur']   = [(in0-1,jn0-1),(in0,jn0-1),(in0,jn0),(in0-1,jn0)]
    sq['ul']   = [(in0,jn0-1),(in0+1,jn0-1),(in0+1,jn0),(in0,jn0)]
    for corner in sq:
        path = sq[corner]
        llpath = Path(toLatlon(latlons,path))
        print 'Corner: ',corner,path,'\n',llpath

        if llpath.contains_point((xlat,xlon)):
            print 'Succes! ',corner,path,llpath
            return path
    print 10*'*','No path found!',xlat,xlon

    if flats[res2] < line_minlat:
        if flons[res2] < line_minlon:
            print 'linksonder'
            i01,j01 = i0,j0
            i02,j02 = i0+1,j0
            i03,j03 = i0+1,j0+1
            i04,j04 = i0,j0+1
        else:
            print 'rechtsonder'
            i01,j01 = i0-1,j0
            i02,j02 = i0,j0
            i03,j03 = i0,j0+1
            i04,j04 = i0-1,j0+1
    else:
        if flons[res2] < line_minlon:
            print 'linksboven'
            i01,j01 = i0-1,j0
            i02,j02 = i0,j0
            i03,j03 = i0,j0+1
            i04,j04 = i0-1,j0+1
        else:
            print 'rechtsboven'
            i01,j01 = i0-1,j0-1
            i02,j02 = i0,j0-1
            i03,j03 = i0,j0
            i04,j04 = i0-1,j0

    altindices = [(j01,i01),(j02,i02),(j03,i03),(j04,i04)]
    print 'zelf: ',altindices,min(altindices)
    return altindices

def minimal_rect_around_line():

    # Create a minimal rectangle of points around the line:
    line_minlat,line_maxlat = min(myslice.line[0]), max(myslice.line[0])
    line_minlon,line_maxlon = min(myslice.line[1]), max(myslice.line[1])

    print line_minlon,line_minlat
    print line_maxlon,line_maxlat
    indicesmin = min(interpol.get_closest_index(myslice.latlons,(line_minlat,line_minlon)))
    indicesmax = max(interpol.get_closest_index(myslice.latlons,(line_maxlat,line_maxlon)))
    print indicesmin,indicesmax

    lat1,lon1 = myslice.line[0][0],myslice.line[1][0]
    lat2,lon2 = myslice.line[0][-1],myslice.line[1][-1]
    altmin = find_four_neighbours( myslice.latlons,(lat1,lon1)) #line_minlat,line_minlon))
    altmax = find_four_neighbours(myslice.latlons,(lat2,lon2)) #(line_maxlat,line_maxlon))
    print 'ff4n min ',altmin
    print 'ff4n max ',altmax

    print 'min',indicesmin,myslice.latlons[0][indicesmin[0],indicesmin[1]],myslice.latlons[1][indicesmin[0],indicesmin[1]]
    print 'max',indicesmax,myslice.latlons[0][indicesmax[0],indicesmax[1]],myslice.latlons[1][indicesmax[0],indicesmax[1]]
    i1,i2 = indicesmin[0],indicesmin[1]
    j1,j2 = indicesmax[0],indicesmax[1]
    print 'i1,i2,j1,j2',i1,i2,j1,j2
    myslice.minrect = [i1,i2,j1,j2]

    print res2, myslice.latlons[0].ravel()[res2], myslice.latlons[1].ravel()[res2]
    altmin.extend(altmax)
    i1,i2 = min([p[0] for p in altmin]), min([p[1] for p in altmin])
    j1,j2 = max([p[0] for p in altmin]), max([p[1] for p in altmin])
    myslice.minrect = [i1,i2,j1,j2] 

    return altmin,altmax


    print 'minlat: ',i1,i2,myslice.latlons[0][i1,i2]
    print 'minlon: ',i1,i2,myslice.latlons[1][i1,i2]
    print 'maxlat: ',j1,j2,myslice.latlons[0][j1,j2]
    print 'maxlon: ',j1,j2,myslice.latlons[1][j1,j2]

    print 'onderste lijn? lon=',i2,myslice.latlons[0][i1:j1,i2],myslice.latlons[1][i1:j1,i2]
    print 'bovenste lijn? lon=',j2,myslice.latlons[0][i1:j1,j2],myslice.latlons[1][i1:j1,j2]

    print 'onderste lijn? lat=',i1,myslice.latlons[0][i1,i2:j2],myslice.latlons[1][i1,i2:j2+1]
    print 'bovenste lijn? lat=',j1,myslice.latlons[0][j1,i2:j2],myslice.latlons[1][j1,i2:j2+1]
    #myslice.minrect = [i1,i2,j1,j2]

def plotmap(domain = 'nl',colors=None,colormap='jet',ml = 60,lsmask_colour = 'black',outdir = './'):

    import matplotlib.pyplot as plt
    import matplotlib.cm as CM
    
    fig  = matplotlib.pyplot.figure()
    ax   = fig.add_subplot(1,1,1)
    
    modelname = 'Harmonie'
    parametername = 'Parameter '+str(myslice.param)
    cmap = CM.get_cmap(colormap)

    plot_map = bmap.myMap(domain = domain, modelname = modelname)
    plot_map.xymap(latlons = (myslice.latlons[0],myslice.latlons[1]), domain = domain, modelname = modelname)
    plot_map.set_axes(ax)
    plot_map.dress_up_map(domain = domain,lsmask_colour = lsmask_colour)
    pcont = plot_map.bmap.contourf(plot_map.x,plot_map.y,myslice.data[ml],# here takes the last, maybe default to level?
                                   20,cmap=cmap,extend='both',ax=ax)
    # colorbar:
    cb = fig.colorbar(pcont,shrink=0.9, extend='both',format='%.2f') #format='%.1e') 
    cb.set_label(parametername,fontsize=9)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(9)

    ax.set_title('Parameter '+str(myslice.param)+' (level '+str(ml)+'), * start, p end', fontsize=9)
    xa,ya   = plot_map.bmap(4.878867,52.375345) #NL148
    print xa,ya
    xs,ys = plot_map.bmap(myslice.line[1],myslice.line[0]) #NL148
    print xs,ys
    plot_map.bmap.plot(xs,ys,'o',color='white')
    plot_map.bmap.plot(xs[0],ys[0],'*',color='white',markersize=15)
    plot_map.bmap.plot(xs[-1],ys[-1],'p',color='white',markersize=15)
        
    plotfile = os.path.join(outdir,'slice_map_'+myslice.case+'_'+str(myslice.param)+'.png')
    fig.savefig(plotfile, dpi = 100)
    print 'Created ',plotfile

if __name__ == '__main__':

    ll1 = (52.3,5.3)
    ll2 = (52.5,3.8)
    myslice.sline(ll1,ll2,30)
    altmin,altmax = minimal_rect_around_line()

    pts = []
    xnpts,ynpts = [],[]
    for p in altmin:
        print p,myslice.latlons[0][p[0],p[1]],myslice.latlons[1][p[0],p[1]]
        pts.append((myslice.latlons[0][p],myslice.latlons[1][p[0],p[1]]))
        xnpts.append(myslice.latlons[0][p])
        ynpts.append(myslice.latlons[1][p])
    scatter(xnpts,ynpts,color='b')

    pts = []
    xpts,ypts = [],[]
    for p in altmax:
        print p,myslice.latlons[0][p[0],p[1]],myslice.latlons[1][p[0],p[1]]
        pts.append((myslice.latlons[0][p],myslice.latlons[1][p[0],p[1]]))
        xpts.append(myslice.latlons[0][p])
        ypts.append(myslice.latlons[1][p])
    scatter(xpts,ypts)

    scatter(ll1[0],ll1[1],color='r')
    scatter(ll2[0],ll2[1],color='r')
    [i1,i2,j1,j2] = myslice.minrect
    xrect = [myslice.latlons[0][i1,i2],myslice.latlons[0][j1,i2],myslice.latlons[0][j1,j2],myslice.latlons[0][i1,j2]]
    yrect = [myslice.latlons[1][i1,i2],myslice.latlons[1][j1,i2],myslice.latlons[1][j1,j2],myslice.latlons[1][i1,j2]]
    scatter(xrect,yrect,color='g')

    #plotmap()
