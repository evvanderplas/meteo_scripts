
if 1:
    # should come from file

    searchlist = []
    levellist  = []
    colorlist  = []
    f = open(settingsfile,'r')
    for line in f:
        pstr = line.split(',')
        print pstr
        a = ('#','\n')
        if any(x in a for x in pstr[0]): 
            print 'dus niet: ',pstr
            pass
        else:
            searchlist.append([pstr[0:3]])

            print pstr[3]
            if isinstance(pstr[3],str):
                # a colormap
                levellist.append([pstr[3]])
            elif isinstance(pstr[3],list):
                # a list of levels
                levellist.append([pstr[3]])
            elif isinstance(pstr[3],tuple):
                try:
                    levellist.append(numpy.arange(pstr[3][0],pstr[3][1],pstr[3][2]))
                except:
                    # only start,stop given, no step, create 10 increments:
                    step = 0.1 * (pstr[3][1] - pstr[3][1])
                    levellist.append(numpy.arange(pstr[3][0],pstr[3][1],step))

            print pstr[4]
            if isinstance(pstr[4],str):
                # a colormap
                colorlist.append(CM.get_cmap(pstr[4]) )
                


    print 'reading from file: '
    print searchlist, levellist, colorlist; #sys.exit(1)
