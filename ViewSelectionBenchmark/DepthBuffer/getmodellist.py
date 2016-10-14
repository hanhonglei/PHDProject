#!/usr/bin/env python
import os, sys 

# list all models in models directory

outfilename = "mlist"
ofile = open(outfilename, 'w') # open file for writing

for f in os.listdir("../humantest/models"): 
    if f != '.svn':
        ofile.write('../humantest/models/%s\n' % f)

ofile.close()
