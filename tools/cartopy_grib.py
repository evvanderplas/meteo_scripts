import os,sys
import cartopy
import matplotlib.pyplot as plt


def main(gfile,p = 'carre'):

    if p == 'carre':
        proj = cartopy.crs.PlateCarree()
    elif p == 'lambert':
        proj = cartopy.crs.LambertConformal(central_longitude=.0, 
                                            central_latitude=60.0, 
                                            false_easting=0.0, 
                                            false_northing=0.0, 
                                            #secant_latitudes=(33, 45), 
                                            secant_latitudes=(50, 60), 
                                            globe=None, cutoff=-30)
    elif p == 'merc':
        proj = cartopy.crs.Mercator(central_longitude=0.0, 
                                    min_latitude=-80.0, 
                                    max_latitude=84.0, 
                                    globe=None)
    elif p == 'eur':
        proj = cartopy.crs.EuroPP()
    #proj = cartopy.crs.PlateCarree()

    ax = plt.axes(projection=proj)
    #ax.stock_img()
    #ax.set_extent([-10, 20, 40, 60]) # lon1,lon2, lat1,lat2
    ax.set_extent([0, 15, 45, 56]) # lon1,lon2, lat1,lat2


    import pygrib
    from matplotlib import cm as CM
    cmap = CM.get_cmap('hot_r')

    
    grbs = pygrib.open(gfile)
    lats,lons = grbs[1].latlons()
    vals      = grbs[1].values
    plt.contourf(lons, lats, vals, [0,.1,.3,1.,5,10,15,25],
                 transform=proj,
                 cmap = cmap)
    plt.colorbar()

    ax.coastlines(resolution='10m')
    cartopy.feature.BORDERS.scale = '10m'
    cartopy.feature.LAKES.scale = '10m'
    ax.add_feature(cartopy.feature.LAKES, linestyle='-')
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
    ax.gridlines()
    
    plt.show()

    
    return 0


    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE, )
    #ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)

    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='gray')
    
    #ax.set_extent([-20, 60, -40, 40])
    ax.set_extent([0, 10, 50, 60]) # lon1,lon2, lat1,lat2
    ax.stock_img()
    
    plt.show()


if __name__ == '__main__':

    gdir = '/usr/people/plas/python/tools' 
    gfile = 'harmonie_400x400.grb'
    gribfile = os.path.join(gdir,gfile)
    main(gribfile,p='eur')
