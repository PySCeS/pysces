#!/usr/bin/env python
# Testing the new parallel scanner class - detection of invalid states
# --jr20101218

import os
backupdir = os.getcwd()

import numpy as np
import pysces

tbox=pysces.PyscesUtils.TimerBox()
import time

m=pysces.model('isola2a')

ser = pysces.Scanner(m)

print "Serial execution..."
print "Start: ", tbox.normal_timer('SER')
print tbox.SER.next()
t1=time.time()
ser.quietRun = True
ser.addScanParameter('V4',0.01,200,2000,log=True)
ser.addUserOutput('J_R1', 'A')
ser.Run()
print "Done: ", tbox.SER.next()
t2=time.time()
print "Duration: %.2f seconds" % (t2-t1)
ser.statespersecond = len(ser.ScanSpace)/(t2-t1)
print "States per second: %.1f" % ser.statespersecond

print "\n\nParallel execution...scans per run =", 36
par = pysces.ParScanner(m)
par.scans_per_run = 36
t3=time.time()
par.addScanParameter('V4',0.01,200,2000,log=True)
par.addUserOutput('J_R1', 'A')
par.Run()
t4=time.time()
print "Duration: %.2f seconds" % (t4-t3)
par.statespersecond = par.Tsteps/(t4-t3)
print "States per second: %.1f" % par.statespersecond

print "\n Speedup: %.2f" % (par.statespersecond/ser.statespersecond)

print "\n\nParallel execution...with scatter and gather"
par2 = pysces.ParScanner(m)
t5=time.time()
par2.addScanParameter('V4',0.01,200,2000,log=True)
par2.addUserOutput('J_R1', 'A')
par2.RunScatter()
t6=time.time()
print "Duration: %.2f seconds" % (t6-t5)
par2.statespersecond = par2.Tsteps/(t6-t5)
print "States per second: %.1f" % par2.statespersecond

print "\n Speedup: %.2f" % (par2.statespersecond/ser.statespersecond)

print "\n===========\nComparing results..."
sss = ser.SteadyStateResults
pss = par.SteadyStateResults
suo = ser.UserOutputResults
puo = par.UserOutputResults
p2ss = par2.SteadyStateResults
p2uo = par2.UserOutputResults

print "\nTesting nan's:"
print "serial and parallel (tc):      ", \
np.alltrue(np.where(np.isnan(sss))[0] == np.where(np.isnan(pss))[0]), \
np.alltrue(np.where(np.isnan(sss))[1] == np.where(np.isnan(pss))[1]), \
np.alltrue(np.where(np.isnan(suo))[0] == np.where(np.isnan(puo))[0]), \
np.alltrue(np.where(np.isnan(suo))[1] == np.where(np.isnan(puo))[1])
print "serial and parallel (scatter): ", \
np.alltrue(np.where(np.isnan(sss))[0] == np.where(np.isnan(p2ss))[0]), \
np.alltrue(np.where(np.isnan(sss))[1] == np.where(np.isnan(p2ss))[1]), \
np.alltrue(np.where(np.isnan(suo))[0] == np.where(np.isnan(p2uo))[0]), \
np.alltrue(np.where(np.isnan(suo))[1] == np.where(np.isnan(p2uo))[1])

print "\nTesting finite values:"
print "serial and parallel (tc):      ", \
np.alltrue(sss[np.where(np.isfinite(sss))]==pss[np.where(np.isfinite(pss))]), \
np.alltrue(suo[np.where(np.isfinite(suo))]==puo[np.where(np.isfinite(puo))])
print "serial and parallel (scatter): ", \
np.alltrue(sss[np.where(np.isfinite(sss))]==p2ss[np.where(np.isfinite(p2ss))]), \
np.alltrue(suo[np.where(np.isfinite(suo))]==p2uo[np.where(np.isfinite(p2uo))])


os.chdir(backupdir)