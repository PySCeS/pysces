"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""

from Kraken import *
from KrakenDataTools import IntersectionAnalysis
os.chdir('c:/mypysces/btk')
##  os.chdir('/home/bgoli/mypysces/kraken')

START = time.time()
"""
Build controller and initialise active tentacles
Get a name for this task (task == collection of jobs)
"""
kcntrl = KrakenScanController()
kcntrl.setEnv(task_id='SCAN_A_P_A0512_big', working_dir='c:\\mypysces\\btk',\
        model_file='isola_fixA.psc', model_dir='c:\\mypysces\\pscmodels\\btk')
kcntrl.startController()
"""
Set up a list of 'static' initialisation commands that will be used to
initialise each server before a job is sent to it
"""
I_list = ['P_INIT,' + kcntrl.kc_model_server.model_file_name, 'P_LOAD',\
'P_SET_FLOAT,nleq2_nitmax,30',\
'P_SET_FLOAT,A05,1.2',\
'P_SET_FLOAT,V2,100.0', 'P_SET_FLOAT,V3,100.0']
kcntrl.setInitCmds(I_list)

pointsP1 = 14
pointsP2 = 1400
scan_segments = 512
##  pointsP1 = 20
##  pointsP2 = 100
##  scan_segments = 20
scan_job = 'P_SCAN, 2, 4, P, %s, %s, '+str(pointsP1)+', 1, 0, A, 0.01, 400,  '+str(pointsP2)+', 1, 0,  P, A, J_R1, J_R3'
kcntrl.setScanJobs(0.01, 40000, scan_segments, scan_job, log=True)

kcntrl.Run(raw_data_dump=False)
kcntrl.concatenateResults()
print('kcntrl.result_array.shape', kcntrl.result_array.shape)

MID = time.time()

kisect = IntersectionAnalysis()

kisect.zmin = 0.01
kisect.setCoords(kisect.DataFilter(kcntrl.getAndClearResultArray(), report=False))
##  kisect.setVTKCoords()

kisect.intersect_tol = 0.04
kisect.CalculateIntersect()
kisect.writeGnuPlot_3D(kisect.intersection, kcntrl.task_id+'_isect')
##  kcntrl.Dump(kisect.intersection, kcntrl.task_id+'_isect.bin')

##  print 'kisect.vtkcoords.shape, type', kisect.vtkcoords.shape, kisect.vtkcoords.dtype

##  kisect.LogVTKdata()
kisect.LogVTKIntersectiondata()
##  kcntrl.Dump(kisect.vtkcoords, kcntrl.task_id+'_vtkcoords.bin')

##  kisect.writeVTK_UnstructuredGrid(kisect.vtkcoords.take([0,1,2],axis=1), kcntrl.task_id+'_J1')
##  kisect.writeVTK_UnstructuredGrid(kisect.vtkcoords.take([0,1,3],axis=1), kcntrl.task_id+'_J3')
kisect.writeVTK_UnstructuredGrid(kisect.vtkintersection.take([0,1,2],axis=1), kcntrl.task_id+'_isect')

END = time.time()

print("\n********** start Kraken report **********")
print("Using %s tentacles" % len(kcntrl.kc_available_server_list))
print("Scan points          = %s" % (pointsP1*pointsP2*scan_segments))
print("Data generation time = %2.2f minutes." % ((MID-START)/60.0))
print("Data analysis time   = %2.2f minutes." % ((END-MID)/60.0))
print("Total scan time      = %2.2f minutes." % ((END-START)/60.0))
print("Total points/minute  = %2.1f" % ((60.0*pointsP1*pointsP2*scan_segments)/(END-START)))
print("********** end Kraken report **********\n")

kcntrl.Dump(kisect, kcntrl.task_id+'_kisect.bin')