import os

humanresult = open("humanresult",'w')

for i in range(0,45) :
    fn = "rmodel-" + str(i)
    f = open(fn,'r')
    vp = ""
    for j in range(0,10) :
        vp = f.readline()
    vp = f.readline()
    humanresult.write(vp)

humanresult.close()



