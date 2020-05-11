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
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import numpy, scipy, itertools, time
try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

from .KrakenNET import socket, time, cPickle
from .KrakenNET import SimpleClient, SimpleMultiReadClient, ThreadedClient
from .KrakenNET import ModelFileServer, TentacleScanner, ServerListLoader
from .KrakenNET import ServerStatusCheck, HOSTNAME, BLOCK_SIZE, PICKLE_PROTOCOL, STATUS_PORT, PYSCES_PORT, CFSERVE_PORT

print('I AM:', HOSTNAME)
print('Using BLOCK_SIZE: %s\n' % BLOCK_SIZE)

class KrakenController:
    kc_status_port = None
    kc_pysces_port = None
    kc_block_size = None

    kc_file_name = 'server_list'
    kc_directory_name = None
    kc_available_server_list = None
    kc_busy_server_list = None
    kc_tentacle_response = None
    kc_model_server = None
    kc_tentacle_scanner = None
    kc_tentacle_monitor = None
    kc_tentacle_process_map = None

    def __init__(self,status_port=STATUS_PORT, pysces_port=PYSCES_PORT, block_size=BLOCK_SIZE):
        self.kc_status_port = status_port
        self.kc_pysces_port = pysces_port
        self.kc_block_size = block_size
        self.kc_available_server_list = []
        self.kc_busy_server_list = []
        self.kc_tentacle_response = {}
        self.kc_tentacle_process_map = {}

        print(self.kc_status_port, self.kc_pysces_port, self.kc_block_size)

    def initialiseTentacles(self,file_name=None, directory_name=os.getcwd()):
        if file_name != None:
            self.kc_file_name = file_name
        self.kc_available_server_list = ServerListLoader().ReadFile(self.kc_file_name, directory_name)
        self.kc_tentacle_scanner = TentacleScanner(self.kc_available_server_list)
        print('Server list:')
        print(self.kc_available_server_list, '\n')
        assert len(self.kc_available_server_list) > 0, '\n* No servers in server list file'

    def getActiveTentacles(self):
        self.kc_tentacle_scanner.scan()
        self.kc_available_server_list = self.kc_tentacle_scanner.servers_ready
        self.kc_busy_server_list = self.kc_tentacle_scanner.servers_busy
        print('Server list:')
        print(self.kc_available_server_list, '\n')
        assert len(self.kc_available_server_list) > 0, '\n* No active servers'
        for t in self.kc_available_server_list:
            self.kc_tentacle_response.setdefault(t)
            self.kc_tentacle_process_map.setdefault(t)
            # maybe wil be needed later for reruns of getActiveTentacles
            ##  if self.kc_tentacle_process_map.has_key[t]:
                ##  dt = self.kc_tentacle_process_map.pop[t]
                ##  del dt

        ##  print self.kc_tentacle_response
        ##  print self.kc_tentacle_process_map


    def Run(self):
        pass

    def sendCmdListToAll(self, command_list, pause=0.5):
        assert type(command_list) == list or type(command_list) == tuple, '\n* Input must be a list or tuple'
        for server in self.kc_available_server_list:
            self.sendCmdListToOne(command_list, server)
            time.sleep(pause)

    def sendCmdListToOne(self, command_list, server):
        assert type(command_list) == list or type(command_list) == tuple, '\n* Input must be a list or tuple'
        assert server in self.kc_available_server_list, '\n* Server %s not in active server list' % server
        clnt = SimpleClient(server, self.kc_pysces_port, self.kc_block_size)
        clnt.send(command_list)
        self.kc_tentacle_response[server] = clnt.response
        ##  print self.kc_tentacle_response

    def sendJobToOne(self, command, server):
        assert server in self.kc_available_server_list, '\n* Server %s not in active server list' % server
        self.kc_tentacle_process_map[server] = ThreadedClient(command, server, self.kc_pysces_port, self.kc_block_size)

    def sendJobToAll(self, command):
        for server in self.kc_available_server_list:
            self.sendJobToOne(command, server)

    def runAllJobs(self):
        Istarted = []
        for s in list(self.kc_tentacle_process_map.keys()):
            if self.kc_tentacle_process_map[s] != None:
                if not self.kc_tentacle_process_map[s].isAlive():
                    self.kc_tentacle_process_map[s].start()
                    Istarted.append(self.kc_tentacle_process_map[s])
        for s in Istarted:
            s.join()
        del Istarted

    def clearProcessMap(self):
        for s in list(self.kc_tentacle_process_map.keys()):
            if self.kc_tentacle_process_map[s] != None:
                if not self.kc_tentacle_process_map[s].isAlive():
                    self.kc_tentacle_process_map[s] = None
        ##  print self.kc_tentacle_process_map


    def getDataFromOneJob(self, server):
        assert server in self.kc_tentacle_process_map, '\n* Server %s not in active server list' % server
        if self.kc_tentacle_process_map[server].response == 'True':
            print('* \"%s\" reports completing job without error' % server)
        else:
            print('* Check result! \"%s\" reports job processing failure'  % server)
        client = SimpleMultiReadClient(server, self.kc_pysces_port, self.kc_block_size)
        client.send('P_GETDATA')
        self.kc_tentacle_response[server] = cPickle.loads(client.response)
        print('\n', server, 'returns', type(self.kc_tentacle_response[server]))
        print('Data successfully retrieved from: %s \n' % server)

    def getDataFromAllJobs(self, pause=0.2):
        for s in list(self.kc_tentacle_process_map.keys()):
            if self.kc_tentacle_process_map[s] != None:
                self.getDataFromOneJob(s)
                time.sleep(pause)

    def killTentacles(self):
        GO = True
        while GO:
            inp = input('\nYou have initiated tentacle shutdown ... are you sure? (y/n): ')
            if inp in ['y','Y','yes']:
                GO = False
                for server in self.kc_available_server_list:
                    SimpleClient(server, self.kc_status_port, self.kc_block_size).send(['KILL'])
                    SimpleClient(server, self.kc_pysces_port, self.kc_block_size).send(['KILL'])
                print('Warning: all servers deleted for this host')
                self.kc_available_server_list = []
            elif inp in ['n','N','no']:
                GO = False
            else:
                print('\nPlease enter y for yes or n for no')

    def chkpsc(self,F):
        try:
            if F[-4:] == '.psc':
                pass
            else:
                print('Assuming extension is .psc')
                F += '.psc'
        except:
            print('Chkpsc error')
        return F

    def startModelServer(self, model_file, model_dir):
        self.model_file = model_file
        self.model_dir = model_dir
        self.model_file = self.chkpsc(self.model_file)
        mpath = os.path.join(self.model_dir, self.model_file)
        if os.path.exists(mpath):
            try:
                self.kc_model_server = ModelFileServer(CFSERVE_PORT, BLOCK_SIZE, 'model_server')
                self.kc_model_server.ReadFile(self.model_file, self.model_dir)
                self.kc_model_server.start()
            except Exception as ex:
                print('Server might already be running on this machine')
                print(ex)
        else:
            print('Error: file %s does not exist!\n' % mpath)
            raise NameError

    def startTentacleMonitor(self, interval=60):
        self.kc_tentacle_monitor = ServerStatusCheck(self.kc_available_server_list, self.kc_status_port, self.kc_block_size, interval)
        self.kc_tentacle_monitor.start()

    def ResetServersToReady(self):
        for server in self.kc_available_server_list:
            SimpleClient(server, self.kc_status_port, self.kc_block_size).send(['P_RESET_STATUS'])


class KrakenScanController(KrakenController):
    job_list = None
    JobsWaiting = None
    result_list = None
    result_array = None
    task_id = None
    working_dir = None
    Model_File = None
    Model_Dir = None
    tentacle_monitor_interval = 120
    init_commands = None
    once_off_tentacle_init = None

    def setEnv(self, task_id, working_dir, model_file, model_dir):
        self.task_id = task_id
        self.working_dir = working_dir
        self.job_list = []
        self.result_list = []
        self.Model_File = model_file
        self.Model_Dir = model_dir

    def startController(self):
        ##  self.blob = getBasicController()
        self.initialiseTentacles(directory_name=self.working_dir)
        self.getActiveTentacles()
        self.startModelServer(self.Model_File, self.Model_Dir)
        self.startTentacleMonitor(interval=self.tentacle_monitor_interval)
        print('Server list: %s\n' % self.kc_available_server_list)

    def setInitCmds(self, lst):
        assert type(lst) == list, 'This must be a list'
        self.init_commands = lst

    def setOnceOffTentacleInit(self, lst):
        assert type(lst) == list, 'This must be a list'
        self.once_off_tentacle_init = lst

    def setJobList(self, lst):
        assert type(lst) == list, 'This must be a list'
        self.job_list = lst


    def setScanJobs(self, start, end, intervals, job, log=False):
        """Splits a range into a number of jobs with intervals"""
        assert intervals >= 1, '\n* Minimum of 1 interval'
        if log:
            kpoints = scipy.logspace(scipy.log10(start), scipy.log10(end), intervals+1)
        else:
            kpoints = scipy.linspace(start, end, intervals+1)
        self.job_list = []
        for p in range(len(kpoints)-1):
            job2 = job % (kpoints[p], kpoints[p+1])
            self.job_list.append(job2)
            print(job2)

    def Run(self, raw_data_dump=True):
        START = time.time()
        self.Scheduler3(self.job_list, self.init_commands)
        MID = time.time()
        if raw_data_dump:
            try:
                self.Dump(self.result_list, self.task_id+'_raw_data')
            except Exception as ex:
                print('Raw data not saved', ex)
        print("Data generation time = %2.2f minutes." % ((MID-START)/60.0))

    def buildCycler(self, lst):
        """
        creates an infinite looping generator with itertools
        """
        return itertools.cycle(lst)

    def Dump(self, thing, fname):
        try:
            F = open(os.path.join(self.working_dir, fname.replace('.bin','')+'.bin'), 'wb')
            cPickle.dump(thing, F, 2)
            F.flush()
            F.close()
        except Exception as ex:
            print('Dump exception raised')
            print(ex)

    def deltaCommand(self):
        """
        This is used to setup conditions for each next-job-on-server
        """
        delta = 'P_NONE'
        return [delta]

    def Scheduler3(self, job_list, init_list):
        """
        job_list is a list of jobs to run
        init_list set of `static' command used to initialise each job
        """
        original_job_list = tuple(job_list)
        TIME_START = time.time()
        self.result_list = []
        self.JobsWaiting = True
        getServer = self.buildCycler(self.kc_available_server_list)
        if self.once_off_tentacle_init != None:
            self.sendCmdListToAll(self.once_off_tentacle_init, pause=0.5)
        while self.JobsWaiting:
            deferred_jobs = []
            job_servers = []
            for job in range(len(job_list)):
                if job < len(self.kc_available_server_list):
                    print('\nJob %s queued for processing...' % (job+1))
                    server = next(getServer)
                    job_servers.append(server)
                    print("\nProcessing job %s init list ..." % (job+1))
                    self.sendCmdListToOne(init_list, server)
                    print("\nProcessing job %s delta command ..." % (job+1))
                    self.sendCmdListToOne(self.deltaCommand(), server)
                    self.sendJobToOne(job_list[job], server)
                else:
                    deferred_jobs.append(job_list[job])
                    print('Job %s deferred and rescheduled.' % (job+1))
            print("\nProcessing queued jobs ...")
            self.runAllJobs()
            for server in job_servers:
                self.getDataFromOneJob(server)
                arr = self.kc_tentacle_response[server]
                self.result_list.append(arr)
                print(type(arr))
            self.clearProcessMap()
            if len(deferred_jobs) > 0:
                print('\nDeferred:', deferred_jobs)
                job_list = deferred_jobs
            else:
                self.JobsWaiting = False

        TIME_END = time.time()

        print('\nNumber of results = %s' % len(self.result_list))
        total_states = 0
        ##  for x in range(len(self.result_list)):
            ##  try:
                ##  print '\tResult %s has %s rows and %s columns' % (x+1, self.result_list[x].shape[0], self.result_list[x].shape[1])
                ##  total_states += self.result_list[x].shape[0]
            ##  except Exception, ex:
                ##  print type(self.result_list[x][0])
                ##  print '\tResult %s %s: %s' % (x+1, type(self.result_list[x][0]), self.result_list[x][0])
                ##  print original_job_list[x]
                ##  print ex
        ##  print 'Total time taken to complete %s state scan %s minutes' % (total_states, (TIME_END-TIME_START)/60.0)

    #investigate using vstack with type/range checking
    def concatenateResults(self):
        output = None
        for arr in range(len(self.result_list)):
            if type(self.result_list[arr]) == numpy.ndarray:
                if arr == 0:
                    output = self.result_list[arr]
                else:
                    output = scipy.concatenate((output,self.result_list[arr]))
        self.result_array = output

    def getResultArray(self):
        return self.result_array

    def getResultList(self):
        return self.result_list

    def getAndClearResultArray(self):
        print('ReSetting result_array and result_list (returning result_array)')
        tmp = self.result_array
        self.result_list = []
        self.result_array = None
        return tmp

    def getAndClearResultList(self):
        print('ReSetting result_array and result_list (returning result_list)')
        tmp = self.result_list
        self.result_list = []
        self.result_array = None
        return tmp



def getBasicController(lvl=3, directory=os.getcwd()):
    if lvl >= 1: blob = KrakenController()
    if lvl >= 2: blob.initialiseTentacles(directory_name=directory)
    if lvl >= 3: blob.getActiveTentacles()
    return blob

