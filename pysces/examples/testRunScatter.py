#!/usr/bin/env python
# Testing the new parallel scanner class

import os
backupdir = os.getcwd()

import numpy as np
import pysces

tbox=pysces.PyscesUtils.TimerBox()
import time

m=pysces.model('isola2a')

print "\n\nParallel execution...using RunScatter"
par2 = pysces.ParScanner(m)
t5=time.time()
par2.addScanParameter('V4',60,100,11)
par2.addScanParameter('V1',100,130,16)
par2.addScanParameter('V2',100,130,16,slave=True)
par2.addScanParameter('V3',80,90,6)
par2.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
#par2.addUserOutput('J_R1', 'A')
par2.RunScatter()
t6=time.time()
print "Duration: %.2f seconds" % (t6-t5)
par2.statespersecond = par2.Tsteps/(t6-t5)
print "States per second: %.1f" % par2.statespersecond

os.chdir(backupdir)