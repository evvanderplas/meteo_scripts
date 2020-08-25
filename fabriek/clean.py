#!/usr/bin/env python

import glob, os
# delete auto generated files

files_to_delete = []

# pyc files
files_to_delete.extend(glob.glob('*.pyc'))
files_to_delete.extend(glob.glob('*/*.pyc'))

# png files
files_to_delete.extend(glob.glob('*.png'))

for f in files_to_delete:
    print 'deleting file: ',f
    os.remove(f)
