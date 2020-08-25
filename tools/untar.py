#! /usr/bin/env python
import os,tarfile
import shutil

def locate(pattern, root=os.curdir):
    '''
    Locate all files matching supplied filename pattern in and below
    supplied root directory.
    '''
    import fnmatch

    # in case pattern is just a string:
    pattern = '*'+pattern+'*'

    for path, dirs, files in os.walk(os.path.abspath(root)):
        #print files
        for filename in fnmatch.filter(files, pattern):
            #print filename
            yield os.path.join(path, filename)

def untar(tfile,outdir,pattern='',tempdir = './temp',verb = False):

    def untar_files(members):
        '''adapted from python website'''
        if '*' in pattern.strip('*'):
            r = pattern[:pattern.find('*')]
        else:
            r = pattern.strip('*')
        for tarinfo in members:
            if r in os.path.split(tarinfo.name)[1]:
                yield tarinfo

    filepath,filename = os.path.split(tfile)
    print '... extracting radar files ...', filename,'with pattern',pattern
    ar = tarfile.open(tfile,'r:*')

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if tempdir == None:
        tempdir = os.path.join(outdir,'temp')
    elif os.path.isdir(tempdir):
        pass
    elif not os.path.isdir(tempdir):
        try:
            os.makedirs(tempdir)
        except OSError:
            print 'Dir exists?',tempdir
            #print 'Hm, no dir, not makeable: ',tempdir; sys.exit(1)

    print 'untarring to ',tempdir
    outsafe = tempdir #os.path.join(outdir,'temp') # os.path.join(outdir,suf2(dtf.year),suf2(dtf.month))
    print filename,' to: ',outsafe

    l = ar.extractall(members=untar_files(ar),path=outsafe)
    print 'using ', pattern,' in locate()',outsafe 
    for xfile in locate(pattern, root=outsafe):
        print xfile
        try:
            shutil.move(xfile,outdir)
        except:
            os.rename(xfile,os.path.join(outdir,os.path.split(xfile)[1]))

    ar.close()
    return 'Done'


if __name__ == '__main__':

    testtar = './somefiles.tar.gz'

    print 'testing untar with',testtar
    untar(testtar,'./newout',tempdir='./newtemp',pattern='file')
    os.listdir('./newout')
