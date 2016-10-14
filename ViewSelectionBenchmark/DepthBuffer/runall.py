import os, sys, string

file = open('mlist', 'r')
idx = 0
for line in file :
    line = line.strip()
    cmd = "dp.exe " + line + " " + str(idx)
    idx += 1
    os.system(cmd)
file.close()
os.system("del v*.pgm v*.bz2")
