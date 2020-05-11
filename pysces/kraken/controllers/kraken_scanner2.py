import os
##  os.chdir('c:/mypysces/kraken')
os.chdir('/home/bgoli/mypysces/kraken')
from Kraken import *
import itertools
import scipy.io

#investigate using vstack with type/range checking
def concatenateArrays(array_list):
    output = None
    for arr in range(len(array_list)):
        if type(array_list[arr]) == numpy.ndarray:
            if arr == 0:
                output = array_list[arr]
            else:
                output = scipy.concatenate((output,array_list[arr]))
    return output

def splitScan(start, end, intervals, job, log=False):
    """Splits a range into a number of jobs with intervals"""
    assert intervals >= 1, '\n* Minimum of 1 interval'
    if log:
        kpoints = scipy.logspace(scipy.log10(start), scipy.log10(end), intervals+1)
    else:
        kpoints = scipy.linspace(start, end, intervals+1)
    job_list = []
    for p in range(len(kpoints)-1):
        job2 = job % (kpoints[p], kpoints[p+1])
        job_list.append(job2)
        print(job2)
    return job_list

def buildCycler(lst):
    """
    creates an infinite looping generator with itertools
    """
    return itertools.cycle(lst)

def deltaCommand():
    """
    This is used to setup conditions for each next-job-on-server
    """
    delta = 'P_NONE'
    return [delta]


def Scheduler1(contrl, job_list, init_list, task_id):
    """
    contrl is a getBasicController() instance
    job_list is a list of jobs to run
    init_list set of `static' command used to initialise each job
    """
    original_job_list = tuple(job_list)
    TIME_START = time.time()
    DaResults = []
    JobsWaiting = True
    getServer = buildCycler(contrl.available_server_list)
    while JobsWaiting:
        deferred_jobs = []
        job_servers = []
        for job in range(len(job_list)):
            if job < len(contrl.available_server_list):
                print('\nJob %s queued for processing...' % (job+1))
                server = next(getServer)
                job_servers.append(server)
                print("\nProcessing job %s init list ..." % (job+1))
                contrl.sendCmdListToOne(init_list, server)
                print("\nProcessing job %s delta command ..." % (job+1))
                contrl.sendCmdListToOne(deltaCommand(), server)
                contrl.sendJobToOne(job_list[job], server)
            else:
                deferred_jobs.append(job_list[job])
                print('Job %s deferred and rescheduled.' % (job+1))
        print("\nProcessing queued jobs ...")
        contrl.runAllJobs()
        for server in job_servers:
            contrl.getDataFromOneJob(server)
            arr = contrl.tentacle_response[server]
            DaResults.append(arr)
            print(type(arr))
        contrl.clearProcessMap()
        if len(deferred_jobs) > 0:
            print('\nDeferred:', deferred_jobs)
            job_list = deferred_jobs
        else:
            JobsWaiting = False
        cPickle.dump(DaResults, open(task_id+'_result_list.bin','wb'), 2)

    TIME_END = time.time()

    print('\nNumber of results = %s' % len(DaResults))
    total_states = 0
    for x in range(len(DaResults)):
        try:
            print('\tResult %s has %s rows and %s columns' % (x+1, DaResults[x][0].shape[0], DaResults[x][0].shape[1]))
            total_states += DaResults[x][0].shape[0]
        except:
            print(type(DaResults[x][0]))
            print('\tResult %s %s: %s' % (x+1, type(DaResults[x][0]), DaResults[x][0]))
            print(original_job_list[x])


    print('Total time taken to complete %s state scan %s minutes' % (total_states, (TIME_END-TIME_START)/60.0))
    return DaResults

def GridSortLR(arr):
    """
    Sort irregular unstructured grid data into monotonically rising grid data
    (where possible) Work from left column to right (uses ultra cool NumPy functions :-)
    """
    print('arr.shape:', arr.shape)
    sorted = numpy.rot90(arr)
    print('rotated.shape:', sorted.shape)
    sorted_idx = numpy.lexsort(sorted)
    print('sorted_idx.shape', sorted_idx.shape)
    sorted = sorted.take(sorted_idx, axis=-1)
    sorted = numpy.flipud(sorted)
    arr = numpy.transpose(sorted)
    print('arr.shape:', arr.shape)
    return arr

def ScalarFunc(val):
    """
    Manipulate the scalar (or z-axis) values
    Accepts and returns a double
    """
    return val

def writeGnuPlot_3D(fname, arr):
    scipy.io.write_array(fname+'_result_gplt.dat', arr, separator=' ')

def writeVTK_UnstructuredGrid(fname, arr):
    """
    Write arr as a VTK unstructured grid to file "fname.vtk"
    arr must either be in:
        x,y,z format (in which case z is used for scalar data) or
        z,y,z,s format (s being point scalar value)
    Overwrite the ScalarFunc(val) function to manipulate the scalar (z-axis) values

    Original version http://mayavi.sourceforge.net/mwiki/PointClouds
    Hacked quite a bit by brett 20070221
    """
    assert arr.shape[1] == 3 or arr.shape[1] == 4, '\nneed s or for columns for this'
    if arr.shape[1] == 4:
        HAVE_SCALARS = 1
    else:
        HAVE_SCALARS = 0
        print('No scalar values supplied Z axis values will be used')

    n=arr.shape[0]
    print("n:",n)
    # write data to vtk polydata file
    # write header
    out = open(fname+'.vtk', 'w')
    h1 = "# vtk DataFile Version 2.0\n"
    h1 += "%s\n" % fname
    h1 += "ASCII\n"
    h1 += "DATASET UNSTRUCTURED_GRID\n"
    h1 += "POINTS " + str(n) + " double\n"
    out.write(h1)
    # write xyz data
    for r in range(n):
        #s = '%15.2f %15.2f %15.2f' % (x[i], y[i], z[i])
        out.write(str(arr[r,0])+" "+str(arr[r,1])+" "+str(arr[r,2])+'\n')

    # write cell data
    out.write("CELLS "+ str(n)+ " "+ str(2*n)+'\n')
    for r in range(n):
            #s = '1 %d \n' % (i)
            out.write("1 "+str(r)+"\n")

    # write cell types
    out.write("CELL_TYPES " + str(n)+'\n')
    for r in range(n):
        out.write("1 \n")

    # write z scalar values
    h2 = '\n' + """POINT_DATA """ + str(n) + "\n"
    h3 = "SCALARS %s double 1\n" % fname
    h3 += "LOOKUP_TABLE default\n"
    out.write(h2 + h3)

    for r in range(n):
        if HAVE_SCALARS:
            sc=(ScalarFunc(arr[r,3]))
        else:
            sc=(ScalarFunc(arr[r,2]))
        out.write(str(sc)+ "\n")

    out.write('\n')
    out.close()




# Parameter scan
if __name__ == '__main__':
    import time
    print(time.strftime('%y%m%d-%H%M'))
    START = time.time()

    """
    Build controller and initialise active tentacles
    """
    blob = getBasicController()
    Model_File = '2006_isola2.psc'
    ##  Model_Dir = 'c:\\mypysces\\pscmodels\\isola06'
    Model_Dir = '/home/bgoli/mypysces/pscmodels/isola06'
    blob.startModelServer(Model_File, Model_Dir)
    blob.startTentacleMonitor(interval=120)
    print('Server list: %s\n' % blob.available_server_list)

    """
    Get a name for this task (task == collection of jobs)
    """
    task_name = 'simple_scan'

    """
    Set up a list of 'static' initialisation commands that will be used to
    initialise each server before a job is sent to it
    """
    I_list = ['P_INIT,' + blob.model_server.model_file_name, 'P_LOAD',\
    'P_SET_FLOAT,mode_state_init,3', 'P_SET_FLOAT,nleq2_nitmax,20',\
    'P_SET_FLOAT,mode_state_init3_factor,0.6',\
    'P_SET_FLOAT,V2,200.0', 'P_SET_FLOAT,V3,200.0']

    pointsP1 = 50
    pointsP2 = 500
    scan_segments = 100
    scan_job = 'P_SCAN, 2, 9, P, %s, %s, '+str(pointsP1)+', 1, 0, A05, 0.01, 10000.0,  '+str(pointsP2)+', 1, 0,  A05, P, J_R1, V2, V3, lambda1, lambda2, eigen_sum, eigen_prod'
    J_list = splitScan(0.01, 3.7e4, scan_segments, scan_job, log=True)


    """
    Feed the controller, job list, initiation list and task name to the scheduler
    and off we go, we end up with a list of results
    """
    DaResults = Scheduler1(blob, J_list, I_list, task_name)
    DaResults = concatenateArrays(DaResults)
    print('\nDaResults.shape', DaResults.shape)
    print('DaResults[0]\n', DaResults[0])

    print(time.strftime('%y%m%d-%H%M'))
    MID = time.time()

    coords = numpy.transpose(numpy.array((DaResults[:,0],DaResults[:,1],DaResults[:,2]),'d'))
    print('\ncoords.shape', coords.shape, coords.dtype)
    print('coords[0]\n', coords[0])
    writeGnuPlot_3D(task_name, coords)
    for r in range(coords.shape[0]):
        coords[r,:] = scipy.log10(coords[r,:])
    writeVTK_UnstructuredGrid(task_name, coords)

    print(time.strftime('%y%m%d-%H%M'))
    END = time.time()

    print("\n********** start Kraken report **********")
    print("Using %s tentacles" % len(blob.available_server_list))
    print("Scan points          = %s" % (pointsP1*pointsP2*scan_segments))
    print("Data generation time = %2.2f minutes." % ((MID-START)/60.0))
    print("Data analysis time   = %2.2f minutes." % ((END-MID)/60.0))
    print("Total scan time      = %2.2f minutes." % ((END-START)/60.0))
    print("Total points/minute  = %2.1f" % ((60.0*pointsP1*pointsP2*scan_segments)/(END-START)))
    print("********** end Kraken report **********\n")



# Continuation
'''
if __name__ == '__main__':
    """
    Build controller and initialise active tentacles
    """
    blob = getBasicController()
    Model_File = '2006_isola2.psc'
    Model_Dir = 'c:\\mypysces\\pscmodels\\isola06'
    ##  Model_Dir = '/home/bgoli/mypysces/pscmodels/isola06'
    blob.startModelServer(Model_File, Model_Dir)
    blob.startTentacleMonitor(interval=60)
    print 'Server list: %s\n' % blob.available_server_list

    """
    Get a name for this task (task == collection of jobs)
    """
    task_name = 'cont_eigen'

    """
    Set up a list of 'static' initialisation commands that will be used to
    initialise each server before a job is sent to it
    """
    I_list = ['P_INIT,' + blob.model_server.model_file_name, 'P_LOAD',\
        'P_SET_FLOAT,mode_state_init,3', 'P_SET_FLOAT,nleq2_nitmax,20',\
        'P_SET_FLOAT,mode_state_init3_factor,0.5',\
        'P_SET_FLOAT,V2,200.0', 'P_SET_FLOAT,V3,200.0', 'P_SET_FLOAT,S05,1.0']
        ##  'P_SET_FLOAT,V2,1000.0', 'P_SET_FLOAT,V3,1000.0', 'P_SET_FLOAT,S05,1.0']

    """
    Set up our continuation:
    We want a parameter to be set with every new job/server and make use
    of the deltaCommand() function. Here we use it to set our y-axis (A05) values
    """
    yrange = 96
    xpoints = 64

    # create a list of values
    delta_list = scipy.logspace(scipy.log10(0.01),scipy.log10(10000.0), yrange)
    # initiate an 'looping' generator with it
    delta_gen = buildCycler(delta_list)

    def deltaCommand():
        """
        This is used to setup conditions for each next-job-on-server
        Used in sendCmdListToOne() so may one or more commands returned
        as a list
        """
        delta = 'P_SET_FLOAT,A05,%s' % delta_gen.next()
        # must return a list !!!
        return [delta]

    """
    Set up our list of jobs, in this case a single continuation+eigen with
    density xpoints at each of the delta (A05 values)
    """
    J_list = []
    for d in delta_list:
        ##  J_list.append('P_CONTINUATION,P,0.001,4.0e2,%s,%s' % (xpoints,d))
        J_list.append('P_CONTINUATION_WITH_EIGEN,P,0.001,3.0e4,%s,%s' % (xpoints,d))

    """
    Feed the controller, job list, initiation list and task name to the scheduler
    and off we go, we end up with a list of results
    """
    DaResults = Scheduler1(blob, J_list, I_list, task_name)

    """
    From here on we are processing/formatting our result data
    """
    iRes = []
    sRes = []
    jRes = []
    eRes = []
    for res in DaResults:
        iRes.append(res[0])
        sRes.append(res[1])
        jRes.append(res[2])
        eRes.append(res[3])
    del DaResults

    FinalIArray = concatenateArrays(iRes)
    FinalSArray = concatenateArrays(sRes)
    FinalJArray = concatenateArrays(jRes)
    FinalEArray = concatenateArrays(eRes)
    print FinalIArray.shape, FinalIArray.dtype
    print FinalSArray.shape, FinalSArray.dtype
    print FinalJArray.shape, FinalJArray.dtype
    print FinalEArray.shape, FinalEArray.dtype
    del iRes
    del sRes
    del jRes
    del eRes

    cPickle.dump((FinalIArray,FinalSArray,FinalJArray,FinalEArray), open(task_name+'_result_array.bin','wb'), 2)

    # get rid of negative values in axes and fluxes
    FinalArrayClean = []
    FinalVTK = []
    for r in range(FinalSArray.shape[0]):
        if FinalIArray[r,:].all() > 0.0 and FinalSArray[r,:].all() > 0.0 and FinalJArray[r,:].all() > 0.0:
            FinalArrayClean.append([FinalIArray[r,0], FinalIArray[r,1], FinalJArray[r,1], max(FinalEArray[r,:].real)])
        else:
            FinalIArray[r,:] == 0.0
            FinalSArray[r,:] == 0.0
            FinalJArray[r,:] == 0.0
            FinalEArray[r,:] == 0.0
    FinalArrayClean = numpy.array(FinalArrayClean)


    # write out an xyz gplt file
    coords = numpy.transpose(numpy.array((FinalArrayClean[:,0],FinalArrayClean[:,1],FinalArrayClean[:,2])))
    # get our data in nice grid format for plotting
    coords = GridSortLR(coords)
    writeGnuPlot_3D(task_name, coords)

    def ScalarFunc(val):
        """
        Manipulate the scalar (z-axis) values
        Accepts and returns a double
        """
        if val >= 0.0:
            return 100.0
        else:
            return val


    # prepare vtk data
    vtkcoords = numpy.transpose(numpy.array((FinalArrayClean[:,0],FinalArrayClean[:,1],FinalArrayClean[:,2],FinalArrayClean[:,3]),'d'))
    for row in range(vtkcoords.shape[0]):
        vtkcoords[row,:3] = scipy.log10(vtkcoords[row,:3])

    print vtkcoords[:10,:]
    vtkcoords = GridSortLR(vtkcoords)
    print vtkcoords[:10,:]

    writeVTK_UnstructuredGrid(task_name+'_result_ug', vtkcoords)
'''
