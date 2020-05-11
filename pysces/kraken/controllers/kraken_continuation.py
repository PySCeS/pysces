import os

os.chdir('c:/mypysces/kraken')
##  os.chdir('/home/bgoli/mypysces/kraken')
from pysces.kraken.Kraken import *
## from Kraken import *

# build controller
blob = getBasicController()
Model_File = '2006_isola2.psc'
Model_Dir = 'c:\\mypysces\\pscmodels\\isola06'
##  Model_Dir = '/home/bgoli/mypysces/pscmodels/isola06'
blob.startModelServer(Model_File, Model_Dir)
blob.startTentacleMonitor(interval=60)

# initialise tentacles
print('Server list: %s\n' % blob.available_server_list)
server_list = blob.available_server_list

I_list = ['P_INIT,' + blob.model_server.model_file_name, 'P_LOAD',\
    'P_SET_FLOAT,mode_state_init,3', 'P_SET_FLOAT,nleq2_nitmax,20',\
    'P_SET_FLOAT,mode_state_init3_factor,0.5',\
    'P_SET_FLOAT,V2,200.0', 'P_SET_FLOAT,V3,200.0']


# continuation
yrange = 100
xpoints = 512
delta_list = scipy.logspace(scipy.log10(0.2),scipy.log10(1.0), yrange)
print(delta_list)
delta_gen = buildCycler(delta_list)

def deltaCommand():
    delta = 'P_SET_FLOAT,A05,%s' % next(delta_gen)
    return [delta]

job_list = []
for d in delta_list:
    ##  job_list.append('P_CONTINUATION,P,0.001,4.0e3,%s,%s' % (xpoints,d))
    job_list.append('P_CONTINUATION,P,0.01,3.0e4,%s,%s' % (xpoints,d))


original_job_list = tuple(job_list)

TIME_START = time.time()
DaResults = []
JobsWaiting = True
getServer = buildCycler(blob.available_server_list)
while JobsWaiting:
    deferred_jobs = []
    job_servers = []
    for job in range(len(job_list)):
        if job < len(server_list):
            server = next(getServer)
            job_servers.append(server)
            blob.sendCmdListToOne(I_list, server)
            blob.sendCmdListToOne(deltaCommand(), server)
            blob.sendJobToOne(job_list[job], server)
        else:
            deferred_jobs.append(job_list[job])
            print('Job %s of %s deferred' % (job+1, len(job_list)))
    blob.runAllJobs()
    ##  blob.sendCmdListToAll(['P_STORE_DATA'])
    for server in job_servers:
        blob.getDataFromOneJob(server)
        arr = blob.tentacle_response[server]
        DaResults.append(arr)
        print(type(arr))
    blob.clearProcessMap()
    if len(deferred_jobs) > 0:
        print('\nDeferred:', deferred_jobs)
        job_list = deferred_jobs
    else:
        JobsWaiting = False
    cPickle.dump(DaResults, open('cont_result_list.bin','wb'), 2)


print('\nNumber of results = %s' % len(DaResults))
total_states = 0
for x in range(len(DaResults)):
    try:
        print('\tResult %s has %s rows and %s columns' % (x+1, DaResults[x].shape[0], DaResults[x].shape[1]))
        total_states += DaResults[x].shape[0]
    except:
        print(type(DaResults[x]))
        print('\tResult %s %s: %s' % (x+1, type(DaResults[x]), DaResults[x]))
        print(original_job_list[x])

TIME_END = time.time()

print('Total time taken to complete %s state scan %s minutes' % (total_states, (TIME_END-TIME_START)/60.0))

FinalArray = concatenateArrays(DaResults)
cPickle.dump(FinalArray, open('cont_result_array.bin','wb'), 2)


if __name__ == '__main__':
    import scipy.io

    # get rid of negative values in axes and fluxes
    FinalArrayClean = []
    FinalVTK = []
    for r in range(FinalArray.shape[0]):
        if FinalArray[r,0] <= 0.0 or FinalArray[r,1] <= 0.0 or FinalArray[r,4] <= 0.0:
            FinalArray[r,:] == 0.0
        else:
            FinalArrayClean.append(FinalArray[r,:])

    FinalArrayClean = numpy.array(FinalArrayClean)

    # write out an xyz gplt file
    coords = numpy.transpose(numpy.array((FinalArrayClean[:,0],FinalArrayClean[:,1],FinalArrayClean[:,4])))
    scipy.io.write_array('cont_result_gplt.dat', coords, separator=' ')

    # write out a f%$^&* VTK file

    """
    Orgiginal version http://mayavi.sourceforge.net/mwiki/PointClouds
    Hacked by brett 20070221

    Converts xyz point cloud into something renderable
    Original code by Rod Holland & Prabhu Ramachandran
    Handy to visualise raw ascii column file
    """

    ##  i=0
    ##  x=[]
    ##  y=[]
    ##  z=[]

    ##  # add log(data)
    ##  for r in range(coords.shape[0]):
        ##  x.append(scipy.log10(coords[r,0]))
        ##  y.append(scipy.log10(coords[r,1]))
        ##  z.append(scipy.log10(coords[r,2]))

    n=coords.shape[0]
    print("n:",n)
    # write data to vtk polydata file
    # write header
    out = open('cont_result_ug.vtk', 'w')
    h1 = "# vtk DataFile Version 2.0\n"
    h1 += "loop\n"
    h1 += "ASCII\n"
    h1 += "DATASET UNSTRUCTURED_GRID\n"
    h1 += "POINTS " + str(n) + " double\n"
    out.write(h1)
    # write xyz data
    for r in range(n):
        #s = '%15.2f %15.2f %15.2f' % (x[i], y[i], z[i])
        out.write(str(scipy.log10(coords[r,0]))+" "+str(scipy.log10(coords[r,1]))+" "+str(scipy.log10(coords[r,2]))+'\n')
        ##  out.write(str(scipy.log10(coords[r,2]))+" "+str(scipy.log10(coords[r,1]))+" "+str(scipy.log10(coords[r,0]))+'\n')

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
    h3 = "SCALARS Z_Value double 1\n"
    h3 += "LOOKUP_TABLE default\n"
    out.write(h2+ h3)

    for r in range(n):
        sc=(scipy.log10(coords[r,2]))
        out.write(str(sc)+ "\n")

    out.write('\n')
    out.close()


    #archived shit
    """
    ##  undup = {}
    ##  rounding = 5
    ##  for r in range(coords.shape[0]):
        ##  undup.setdefault(round(coords[r,2], rounding), [coords[r,0], coords[r,1]])
    ##  rF = open('cont_result_array_gplt_reduced.dat','w')
    ##  print 'Total states: %s\nReduced states: %s ' % (total_states, len(undup.keys()))
    ##  for key in undup:
        ##  rF.write('%s, %s, %s\n' % (undup[key][0], undup[key][1], key))
    ##  rF.flush()
    ##  rF.close()

    """
