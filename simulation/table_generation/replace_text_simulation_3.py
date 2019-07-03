#!/usr/bin/env python3
import fileinput

f = open("simulation_3.txt",'r')
filedata = f.read()
f.close()

newdata = filedata
newdata = newdata.replace("1\_bias","bias")
newdata = newdata.replace("1\_cp","cp")
newdata = newdata.replace("1\_sd","sd")
newdata = newdata.replace("2\_bias","bias")
newdata = newdata.replace("2\_cp","cp")
newdata = newdata.replace("2\_sd","sd")
newdata = newdata.replace("3\_bias","bias")
newdata = newdata.replace("3\_cp","cp")
newdata = newdata.replace("3\_sd","sd")
newdata = newdata.replace("4\_bias","bias")
newdata = newdata.replace("4\_cp","cp")
newdata = newdata.replace("4\_sd","sd")

newdata = newdata.replace("alpha=0","$\alpha=0$")

newdata = newdata.replace("n &","$n$ &")
newdata = newdata.replace("dgm\_type","GM")

newdata = newdata.replace("backslash 2","backslash 6")

f = open("simulation_3_modified.txt",'w')
f.write(newdata)
f.close()