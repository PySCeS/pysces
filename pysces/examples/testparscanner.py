#!/usr/bin/env python
# Testing the new parallel scanner class

import os

backupdir = os.getcwd()

import numpy as np
import pysces

tbox = pysces.PyscesUtils.TimerBox()
import time

m = pysces.model('isola2a')

ser = pysces.Scanner(m)

print("Serial execution...")
print("Start: ", tbox.normal_timer('SER'))
print(next(tbox.SER))
t1 = time.time()
ser.quietRun = True
ser.addScanParameter('V4', 60, 100, 11)
ser.addScanParameter('V1', 100, 160, 16)
ser.addScanParameter('V2', 100, 130, 16, follower=True)
ser.addScanParameter('V3', 80, 90, 6)
ser.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# ser.addUserOutput('J_R1', 'A')
ser.Run()
print("Done: ", next(tbox.SER))
t2 = time.time()
print("Duration: %.2f seconds" % (t2 - t1))
ser.statespersecond = len(ser.ScanSpace) / (t2 - t1)
print("States per second: %.1f" % ser.statespersecond)

print("\n\nParallel execution...scans per run =", 100)
par = pysces.ParScanner(m, engine='multiproc')
par.scans_per_run = 100
t3 = time.time()
par.addScanParameter('V4', 60, 100, 11)
par.addScanParameter('V1', 100, 160, 16)
par.addScanParameter('V2', 100, 130, 16, follower=True)
par.addScanParameter('V3', 80, 90, 6)
par.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# par.addUserOutput('J_R1', 'A')
par.Run()
t4 = time.time()
print("Duration: %.2f seconds" % (t4 - t3))
par.statespersecond = par.Tsteps / (t4 - t3)
print("States per second: %.1f" % par.statespersecond)

print(
    "\n Speedup with load balanced TaskClient: %.2f"
    % (par.statespersecond / ser.statespersecond)
)

print("\n\nParallel execution...using RunScatter")
par2 = pysces.ParScanner(m, engine='ipcluster')
t5 = time.time()
par2.addScanParameter('V4', 60, 100, 11)
par2.addScanParameter('V1', 100, 160, 16)
par2.addScanParameter('V2', 100, 130, 16, follower=True)
par2.addScanParameter('V3', 80, 90, 6)
par2.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# par2.addUserOutput('J_R1', 'A')
par2.RunScatter()
t6 = time.time()
print("Duration: %.2f seconds" % (t6 - t5))
par2.statespersecond = par2.Tsteps / (t6 - t5)
print("States per second: %.1f" % par2.statespersecond)

print("\n Speedup with RunScatter: %.2f" % (par2.statespersecond / ser.statespersecond))

print("\n===========\nComparing results...")
# comp = np.equal(ser.SteadyStateResults, par.SteadyStateResults[::5])
# print np.alltrue(comp)
# comp2 = np.equal(ser.UserOutputResults, par.UserOutputResults[::5])
# print np.alltrue(comp)

sss = ser.SteadyStateResults
pss = par.SteadyStateResults
suo = ser.UserOutputResults
puo = par.UserOutputResults
p2ss = par2.SteadyStateResults
p2uo = par2.UserOutputResults

print("serial vs. parallel s/s results     : ", np.alltrue(np.equal(sss, pss)))
print("serial vs. parallel user output     : ", np.alltrue(np.equal(suo, puo)))
print("TaskClient vs RunScatter s/s results: ", np.alltrue(np.equal(pss, p2ss)))
print("TaskClient vs RunScatter user output: ", np.alltrue(np.equal(puo, p2uo)))


os.chdir(backupdir)
