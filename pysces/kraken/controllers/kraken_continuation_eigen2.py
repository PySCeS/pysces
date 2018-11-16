from Kraken import *
from KrakenDataTools import DataTools

os.chdir('c:/mypysces/btk')
##  os.chdir('/home/bgoli/mypysces/kraken')

START = time.time()
"""
Build controller and initialise active tentacles
Get a name for this task (task == collection of jobs)
"""
kcntrl = KrakenScanController()
kcntrl.setEnv(task_id='isola_cont_eigen', working_dir='c:\\mypysces\\btk',\
        model_file='2006_isola1_dem.psc', model_dir='c:\\mypysces\\pscmodels\\isola06')
kcntrl.startController()
"""
Set up a list of 'static' initialisation commands that will be used to
initialise each server before a job is sent to it
"""
I_list = ['P_INIT,' + kcntrl.kc_model_server.model_file_name, 'P_LOAD',\
    'P_SET_FLOAT,mode_state_init,3', 'P_SET_FLOAT,nleq2_nitmax,30',\
    'P_SET_FLOAT,mode_state_init3_factor,0.5', 'P_SET_FLOAT,pitcon_iter,20',\
    'P_SET_FLOAT,V2,100.0', 'P_SET_FLOAT,V3,100.0', 'P_SET_FLOAT,S05,1.0']
    ##  'P_SET_FLOAT,V2,1000.0', 'P_SET_FLOAT,V3,1000.0', 'P_SET_FLOAT,S05,1.0']
kcntrl.setInitCmds(I_list)


"""
Set up our continuation:
We want a parameter to be set with every new job/server and make use
of the deltaCommand() function. Here we use it to set our y-axis (A05) values
"""
# node split on yrange
# A05
yrange = 57
# P start points
xpoints = 200

# create a list of values A05
##  delta_list_vl = scipy.logspace(scipy.log10(0.1),scipy.log10(5.0), yrange)
##  delta_list_l = scipy.logspace(scipy.log10(5.1),scipy.log10(100.0), yrange)
##  delta_list_h = scipy.logspace(scipy.log10(100.1),scipy.log10(5000.0), yrange*3)
delta_list = scipy.logspace(scipy.log10(0.5),scipy.log10(1000.0), yrange)

# initiate an 'looping' generator with it
delta_gen = kcntrl.buildCycler(delta_list)

def deltaCommand():
    """
    This is used to setup conditions for each next-job-on-server
    Used in sendCmdListToOne() so may be one or more commands returned
    as a list
    """
    delta = 'P_SET_FLOAT,A05,%s' % next(delta_gen)
    # must return a list !!!
    return [delta]

kcntrl.deltaCommand = deltaCommand

"""
Set up our list of jobs, in this case a single continuation+eigen with
density xpoints at each of the delta (A05 values)
"""
J_list = []
for d in delta_list:
    ##  J_list.append('P_CONTINUATION,P,0.001,4.0e2,%s,%s' % (xpoints,d))
    J_list.append('P_CONTINUATION_WITH_EIGEN,P,0.01,3.8e4,%s,%s' % (xpoints,d))

kcntrl.setJobList(J_list)
kcntrl.Run(raw_data_dump=False)


DaResults = kcntrl.getAndClearResultList()

# Have data will travel

kdata = DataTools()

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

def concatenateArrays(array_list):
    output = None
    for arr in range(len(array_list)):
        if type(array_list[arr]) == numpy.ndarray:
            if arr == 0:
                output = array_list[arr]
            else:
                output = scipy.concatenate((output,array_list[arr]))
    return output

FinalIArray = concatenateArrays(iRes)
print('Indexes[0]', FinalIArray[0])
FinalSArray = concatenateArrays(sRes)
FinalJArray = concatenateArrays(jRes)
FinalEArray = concatenateArrays(eRes)
print('FinalIArray', FinalIArray.shape, FinalIArray.dtype)
print('FinalSArray', FinalSArray.shape, FinalSArray.dtype)
print('FinalJArray', FinalJArray.shape, FinalJArray.dtype)
print(FinalEArray.shape, FinalEArray.dtype)
del iRes
del sRes
del jRes
del eRes

kcntrl.Dump((FinalIArray,FinalSArray,FinalJArray,FinalEArray), kcntrl.task_id+'_result_array.bin')

# get rid of negative values in axes and fluxes
FinalArrayCleanJ1 = []
FinalArrayCleanJ4 = []
FinalVTK = []
for r in range(FinalSArray.shape[0]):
    if FinalIArray[r,:].all() > 0.0 and FinalSArray[r,:].all() > 0.0 and \
        FinalJArray[r,:].all() > 0.0 and FinalIArray[r,1] >= 0.2 and FinalJArray[r,1] >= 0.08 and FinalIArray[r,1] <= 40000.0:
        FinalArrayCleanJ1.append([FinalIArray[r,0], FinalIArray[r,1], FinalJArray[r,0], max(FinalEArray[r,:].real)])
        FinalArrayCleanJ4.append([FinalIArray[r,0], FinalIArray[r,1], FinalJArray[r,3], max(FinalEArray[r,:].real)])
    else:
        FinalIArray[r,:] = 0.0
        FinalSArray[r,:] = 0.0
        FinalJArray[r,:] = 0.0
        FinalEArray[r,:] = 0.0
FinalArrayCleanJ1 = numpy.array(FinalArrayCleanJ1)
FinalArrayCleanJ4 = numpy.array(FinalArrayCleanJ4)

kdata.writeGnuPlot_3D(FinalArrayCleanJ1, kcntrl.task_id+'_J1')

def ScalarFunc(val):
    """
    Manipulate the scalar (z-axis) values
    Accepts and returns a double
    """
    if val >= 0.0:
        return 60.0
    else:
        return val


kdata.vtk_scalar_func = ScalarFunc

# prepare vtk data
vtkcoords = numpy.transpose(numpy.array((FinalArrayCleanJ1[:,0],FinalArrayCleanJ1[:,1],FinalArrayCleanJ1[:,2],FinalArrayCleanJ1[:,3]),'d'))
for row in range(vtkcoords.shape[0]):
    vtkcoords[row,:3] = scipy.log10(vtkcoords[row,:3])

print(vtkcoords[:5,:])
vtkcoords = kdata.GridSortLR(vtkcoords)
print(vtkcoords[:5,:])

kdata.writeVTK_UnstructuredGrid(vtkcoords, kcntrl.task_id+'_result_ug')
