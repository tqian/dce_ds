#!/usr/bin/env python3
import fileinput

f = open("data_analysis_1.txt",'r')
filedata = f.read()
f.close()

newdata = filedata
newdata = newdata.replace("est1\_1.mortality","mort")
newdata = newdata.replace("est1\_2.ci","ci")
newdata = newdata.replace("est2\_1.mortality","mort")
newdata = newdata.replace("est2\_2.ci","ci")
newdata = newdata.replace("est3\_1.mortality","mort")
newdata = newdata.replace("est3\_2.ci","ci")
newdata = newdata.replace("est4\_1.mortality","mort")
newdata = newdata.replace("est4\_2.ci","ci")

newdata = newdata.replace("1.1061","$\infty$")
newdata = newdata.replace("2.730","2")
newdata = newdata.replace("3.547","1.5")
newdata = newdata.replace("4.365","1")

newdata = newdata.replace("gamma","$\gamma$")

f = open("data_analysis_1_modified.txt",'w')
f.write(newdata)
f.close()