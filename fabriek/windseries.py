#! /usr/bin/env python
import os,sys,glob,pickle,time
import numpy,pygrib
import datetime as dt

def read_station(afile):

    print 'opening',afile
    f = open(afile,'r')

    station = []

    for line in f:
        ks = line.split()
        try: 
            if ks[0] == '#' and ks[1].isdigit:
                if ks[1] == 'STN':
                    print ks[1],ks[2],ks[3],ks[4]
                    station.append([ks[1][0:3],(float(ks[2]),float(ks[3]),float(ks[4])/10)])

            else:
                print 'klaar'
                break
        except IndexError:
            print 'lijst afgelopen'
            break

    f.close()

    return station

def read_ascii_data(afile,ts_object,dtg = '2012062'):

    print 'opening',afile
    f = open(afile,'r')

    # read raw file:
    for line in f:

        ks = line.split()
   
        station = dict()
        # first read station locations from file
        if len(ks) > 1 and ks[0] == '#': 
            print 'preliminaries: ',ks

            if ks[1] == 'STN':

                station = read_station(afile)
                station = dict(station)
                print station

            elif 'RH,' in ks:
                legend = ks            
                if 'RH,' in legend:
                    i_RH = legend.index('RH,')
                else:
                    pass
        
                print 'waar is RH? ',legend,i_RH; #sys.exit(1)

                if 'FF,' in legend:
                    i_F = legend.index('FF,')
                elif 'FX,' in legend:
                    i_F = legend.index('FX,')
                else:
                    print 'no wind ', legend


                if 'T,' in legend:
                    i_T = legend.index('T,')
                elif 'TX,' in legend:
                    i_T = legend.index('TX,')
                else:
                    print 'no temperature', legend

        elif len(ks) == 1: pass #ks[0] != '#': sys.exit(1)

        # read station data
        else:

            # now the data is split by comma's: ',':
            ks = line.split(',')
            #print 'en nu:', ks[1],ks[i_T], ks[i_RH], ks[i_F]; #sys.exit(1)

            if  dtg in ks[1] and ks[7].strip().isdigit(): # and int(ks[7]) > 1:
            
                HH = ks[2] # for hourly stat
                #HH = '0'     # for daily stat

                # message_type: 'ADPSFC'
                st   = int(ks[0])                      # station ID
                dt   = int(ks[1])
                h    = int(ks[2])

                print ks[0].strip()
                print station

                #time = ks[1]+'_'+suf2(HH)+'0000'
                #lat  = float(station[ks[0].strip()][1])
                #lon  = float(station[ks[0].strip()][0])
                #alt  = float(station[ks[0].strip()][2])  # altitude
                # gribcode, acc.interval
                try:
                    RH   = float(ks[i_RH])/10            # value
                    T2M  = float(ks[i_T])/10            # value
                    FF   = float(ks[i_F])/10            # value
                except:
                    print 'Formatting error, passing:'
                    RH   = 0.
                    T2M  = 0.
                    FF   = 0.
                finally:
                    #aug2006.addobsdata(240,11,20060801,00,12.8)
                    print 'writing to file ',st,dt,h,T2M,FF,RH
                    ts_object.addobsdata(st,11,dt,h,T2M)
                    ts_object.addobsdata(st,33,dt,h,FF)
                    ts_object.addobsdata(st,61,dt,h,RH)

    f.close()
    return ts_object


if __name__ == '__main__':

    # setting up parameter list to plot
    import param_lists,settings,plottypes
    import tools.read_stations as read_stations
    import tools.timeseries as ts
    import sources

    stationfile = 'tools/stations.txt'
    obsfile     = 'tools/KNMI_20120701_hourly.txt'
    ts_dict = dict()
    tsfile  = './serieswinddata.pkl' 

    # generate a list of (grib)files for the specified sources
    filelist = []
    for source in sources.sources_list:
        print source['name']
        a = sorted( glob.glob( os.path.join(source['dir'],source['rexp'])) )
        filelist.append( a )

    #print filelist

    t0 = time.clock()

    for s,source in enumerate(sources.sources_list):

        # could be file-source dependent (?)
        outdir = './'
        if source['shortname'] == 'EXP': outdir = source['outdir']
        #print outdir; sys.exit(1)

        stations    = read_stations.stationdict(stationfile,filelist[s][0], model = source)

        tsfile = 'ts_'+source['shortname']+'.pkl'
        if not os.path.isfile(tsfile):
            
            ts_dict[source['shortname']] = ts.stationlist()
            for station in stations:
                stat = stations[station]
                ts_dict[source['shortname']].add_station(stat['name'],
                                                         stat['latlon'],
                                                         stat['stationid'])

            # add observation data to file
            ts_dict[source['shortname']] = read_ascii_data(obsfile,ts_dict[source['shortname']])
            
            print 'write to file ',tsfile
            fp = open(tsfile,'wb')
            pickle.dump(ts_dict[source['shortname']],fp)
            fp.close()

        else:
            print 'opening file ',tsfile
            fp = open(tsfile,'r')
            ts_dict[source['shortname']] = pickle.load(fp)
            fp.close()


        date1 = dt.datetime(2012,6,1,0)
        date2 = dt.datetime(2012,6,22,0)
        print 'get obsdata, ',ts_dict[source['shortname']].stationdatadict[240].get_obsdata(11,20120621,20120622)
        print ts_dict[source['shortname']].stationdatadict[240].get_data(11,date1=20120621,date2=20120622)
        print ts_dict[source['shortname']].stationdatadict[240].obsdict[11].obsdatalist[1].obsdata       # punt: value of date
        sys.exit(1)

        grbindx_prev = None
        for grib_file in filelist[s][:2]:
            
            print '************ open file ',grib_file
            print 'processing ',source['name'],time.clock() - t0, int(time.clock() - t0)*'='

            grbindx=pygrib.index(grib_file,
                                 'indicatorOfParameter',
                                 'indicatorOfTypeOfLevel',
                                 'level')
            
            plottypes.point_extract(grbindx, source, stations, ts_dict[source['shortname']], outdir = outdir)
            #sys.exit(1)

    




















