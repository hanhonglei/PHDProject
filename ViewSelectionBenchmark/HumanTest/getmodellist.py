#!/usr/bin/env python
import os, sys 

# list all models in models directory

outfilename = "modellist"
ofile = open(outfilename, 'w') # open file for writing

for f in os.listdir("./models"): 
    if f != ".svn":
        ofile.write('./models/%s\n' % f)

ofile.close()
