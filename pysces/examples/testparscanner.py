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
ser.addScanParameter('V4', 60, 100, 41)
ser.addScanParameter('V1', 100, 160, 31)
ser.addScanParameter('V2', 100, 130, 31, follower=True)
ser.addScanParameter('V3', 80, 90, 11)
ser.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# ser.addUserOutput('J_R1', 'A')
ser.Run()
print("Done: ", next(tbox.SER))
t2 = time.time()
print("Duration: %.2f seconds" % (t2 - t1))
ser.statespersecond = len(ser.ScanSpace) / (t2 - t1)
print("States per second: %.1f" % ser.statespersecond)

print("\n\nParallel execution...multiproc...scans per run =", 874)
par = pysces.ParScanner(m, engine='multiproc')
par.scans_per_run = 874
t3 = time.time()
par.addScanParameter('V4', 60, 100, 41)
par.addScanParameter('V1', 100, 160, 31)
par.addScanParameter('V2', 100, 130, 31, follower=True)
par.addScanParameter('V3', 80, 90, 11)
par.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# par.addUserOutput('J_R1', 'A')
par.Run()
t4 = time.time()
print("Duration: %.2f seconds" % (t4 - t3))
par.statespersecond = par.Tsteps / (t4 - t3)
print("States per second: %.1f" % par.statespersecond)

print(
    "\n Speedup with load balanced TaskClient (multiproc): %.2f"
    % (par.statespersecond / ser.statespersecond)
)

print("\n\nParallel execution...ipcluster...scans per run =", 874)
par3 = pysces.ParScanner(m, engine='ipcluster')
par3.scans_per_run = 874
t3 = time.time()
par3.addScanParameter('V4', 60, 100, 41)
par3.addScanParameter('V1', 100, 160, 31)
par3.addScanParameter('V2', 100, 130, 31, follower=True)
par3.addScanParameter('V3', 80, 90, 11)
par3.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# par3.addUserOutput('J_R1', 'A')
par3.Run()
t4 = time.time()
print("Duration: %.2f seconds" % (t4 - t3))
par3.statespersecond = par3.Tsteps / (t4 - t3)
print("States per second: %.1f" % par3.statespersecond)

print(
    "\n Speedup with load balanced TaskClient (ipcluster): %.2f"
    % (par3.statespersecond / ser.statespersecond)
)

print("\n\nParallel execution...using RunScatter")
par2 = pysces.ParScanner(m, engine='ipcluster')
t5 = time.time()
par2.addScanParameter('V4', 60, 100, 41)
par2.addScanParameter('V1', 100, 160, 31)
par2.addScanParameter('V2', 100, 130, 31, follower=True)
par2.addScanParameter('V3', 80, 90, 11)
par2.addUserOutput('J_R1', 'A', 'ecR4_X', 'ccJR1_R1')
# par2.addUserOutput('J_R1', 'A')
par2.RunScatter()
t6 = time.time()
print("Duration: %.2f seconds" % (t6 - t5))
par2.statespersecond = par2.Tsteps / (t6 - t5)
print("States per second: %.1f" % par2.statespersecond)

print("\n Speedup with RunScatter: %.2f" % (par2.statespersecond / ser.statespersecond))

print("\n===========\nComparing results...")
sss = ser.SteadyStateResults
pss = par.SteadyStateResults
suo = ser.UserOutputResults
puo = par.UserOutputResults
p2ss = par2.SteadyStateResults
p2uo = par2.UserOutputResults
p3ss = par3.SteadyStateResults
p3uo = par3.UserOutputResults

print("serial vs. parallel s/s results     : ", np.all(np.equal(sss, pss)))
print("serial vs. parallel user output     : ", np.all(np.equal(suo, puo)))
print("TaskClient vs RunScatter s/s results: ", np.all(np.equal(pss, p2ss)))
print("TaskClient vs RunScatter user output: ", np.all(np.equal(puo, p2uo)))
print("multiproc vs ipcluster s/s results  : ", np.all(np.equal(pss, p3ss)))
print("multiproc vs ipcluster user output  : ", np.all(np.equal(puo, p3uo)))


os.chdir(backupdir)
