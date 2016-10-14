#!/usr/bin/env python
import os, sys 

# list all models in models directory

outfilename = "rlist"
ofile = open(outfilename, 'w') # open file for writing

for f in os.listdir("."): 
    if f[0:8] == 'view.dat' :
        ofile.write('%s\n' % f)
ofile.close()
