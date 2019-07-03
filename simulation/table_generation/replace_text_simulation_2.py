#!/usr/bin/env python3
import fileinput

f = open("simulation_2.txt",'r')
filedata = f.read()
f.close()

newdata = filedata
newdata = newdata.replace("1\_1\_bias","bias")
newdata = newdata.replace("1\_1\_cp","cp(\%)")
newdata = newdata.replace("1\_1\_sd","sd")
newdata = newdata.replace("3\_1\_bias","bias")
newdata = newdata.replace("3\_1\_cp","cp(\%)")
newdata = newdata.replace("3\_1\_sd","sd")
newdata = newdata.replace("1\_bias","bias")
newdata = newdata.replace("1\_cp","cp(\%)")
newdata = newdata.replace("1\_sd","sd")
newdata = newdata.replace("3\_bias","bias")
newdata = newdata.replace("3\_cp","cp(\%)")
newdata = newdata.replace("3\_sd","sd")

newdata = newdata.replace("n &","$n$ &")
newdata = newdata.replace("dgm\_type","GM")

newdata = newdata.replace("backslash 2","backslash 6")


f = open("simulation_2_modified.txt",'w')
f.write(newdata)
f.close()