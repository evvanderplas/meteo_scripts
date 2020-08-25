#! /usr/bin/env python

def suf2(hh):
    h = int(hh)
    if   h==0: return '00'
    elif h<10: return '0'+str(h)
    elif h<100: return str(h)
    elif h<1000: return '0'+str(h/100)
    elif h<10000: return str(h/100)
    else: 
        print 'suf2 not well implemented: ',hh,h
        sys.exit(0)

def suf3(h):
    return '0'+suf2(h)

if __name__ == '__main__':

    print 'create prefix zeros: ',3,suf2(3),'13',suf3('13')

    print 'using zfill: ',3,'3'.zfill(2),'3'.zfill(3)
