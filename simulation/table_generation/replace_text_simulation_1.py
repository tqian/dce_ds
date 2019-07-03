#!/usr/bin/env python3
import fileinput

f = open("simulation_1.txt",'r')
filedata = f.read()
f.close()

newdata = filedata
newdata = newdata.replace("12\_bias","bias")
newdata = newdata.replace("12\_cp","cp(\%)")
newdata = newdata.replace("12\_sd","sd")
newdata = newdata.replace("13\_bias","bias")
newdata = newdata.replace("13\_cp","cp(\%)")
newdata = newdata.replace("13\_sd","sd")
newdata = newdata.replace("1\_bias","bias")
newdata = newdata.replace("1\_cp","cp(\%)")
newdata = newdata.replace("1\_sd","sd")
newdata = newdata.replace("3\_bias","bias")
newdata = newdata.replace("3\_cp","cp(\%)")
newdata = newdata.replace("3\_sd","sd")
newdata = newdata.replace("7\_bias","bias")
newdata = newdata.replace("7\_cp","cp(\%)")
newdata = newdata.replace("7\_sd","sd")


newdata = newdata.replace("n &","$n$ &")
newdata = newdata.replace("dgm\_type","GM")

newdata = newdata.replace("backslash 2","backslash 6")

f = open("simulation_1_modified.txt",'w')
f.write(newdata)
f.close()