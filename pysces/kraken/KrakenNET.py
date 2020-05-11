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

import os, time
import pickle, io
import socket
from threading import Thread
from .startup import octopussy

print(octopussy)

HOSTNAME = socket.gethostbyaddr(socket.gethostname())[0]
BLOCK_SIZE = 32768
PICKLE_PROTOCOL = 2
STATUS_PORT = 60000
PYSCES_PORT = 60001
CFSERVE_PORT = 60005

class SimpleClient:
    """
    Sends a list of commands: data = [cmd1,cmd2, ...]
    Standard blocking IO where for each cmd sent a response
    is expected and held in self.response. This is meant
    for quick sequential jobs (housekeeping functions etc)
    """
    server = None
    port = None
    block_size = None
    sent = None
    response = None
    timeout = None

    def __init__(self, server, port, block_size, myname=None):
        self.server = server
        self.port = port
        self.block_size = block_size
        self.sent = []
        self.response = []

    def send(self, data):
        self.sent = data
        self.response = []
        for d in self.sent:
            self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.s.settimeout(self.timeout)
            self.s.connect((self.server, self.port))
            self.s.send(d)
            print('Sent: ', d)
            self.response.append(self.s.recv(self.block_size))
            self.s.close()


class SimpleMultiReadClient:
    """
    Reads a data block which is larger than block_size
    Standard blocking IO which issues a single command (data)
    and reads block_size until all data is returned

    """
    server = None
    port = None
    block_size = None
    sent = None
    response = None
    timeout = None

    def __init__(self, server, port, block_size, myname=None):
        self.server = server
        self.port = port
        self.block_size = block_size

    def send(self, data):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.settimeout(self.timeout)
        self.s.connect((self.server, self.port))

        self.s.send(data)
        print('Sent: ', data)
        GO = True
        self.response = ''
        while GO:
            data = self.s.recv(self.block_size)
            self.response += data
            if data == '':
                GO = False
        self.s.close()


class ThreadedClient(Thread):
    """
    Standard blocking IO where command is a single P_PROTOCOL command which
    expects a response (held in self.response), however, every instance
    of this client is a new thread. This is for long running commands (jobs)
    """
    server = None
    port = None
    block_size = None
    response = None
    timeout = None

    def __init__(self, command, server, port, block_size, myname=None):
        Thread.__init__(self)
        self.command = command
        self.server = server
        self.port = port
        self.block_size = block_size
        if myname != None:
            self.setName(myname)

    def run(self):
        self.SendLoop()

    def SendLoop(self):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.settimeout(self.timeout)
        self.s.connect((self.server, self.port))
        self.s.send(self.command)
        self.response = self.s.recv(self.block_size)
        print(self.getName() + ' Sent:    ', self.command)
        self.s.close()
        #time.sleep(0.5)

class BasicServerSocket(Thread):
    backlog = 5
    port = None
    block_size = None
    client = None
    client_address = None
    RequestLogOn = True
    server_active = True

    def __init__(self, port, block_size, myname=None):
        Thread.__init__(self)
        self.port = port
        self.block_size = block_size
        if myname != None:
            self.setName(myname)
        print(self.getName() + ': Ready to serve!')

    def run(self):
        self.ListenLoop()

    def SendAction(self, data):
        print(self.client_address[0] + ', ' + time.strftime('%H:%M:%S') + ', ' + self.getName() + ', ' + data)
        return data

    def RequestLog(self, data):
        print(self.client_address[0] + ', ' + time.strftime('%H:%M:%S') +\
            ', ' + self.getName() + ', ' + data)

    def KillCheck(self, data):
        if data[:4] == 'KILL':
            data = 'You killed the server at: ' + time.strftime('%H:%M:%S')
            self.server_active = False
            print(self.client_address[0] + ' terminated me (' + self.getName() + ') at '+ time.strftime('%H:%M:%S'))
        return data

    def ListenLoop(self):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.bind(('', self.port))
        self.s.listen(self.backlog)
        while self.server_active:
            self.client, self.client_address = self.s.accept()
            data = self.client.recv(self.block_size)
            if self.RequestLogOn: self.RequestLog(data)
            data = self.KillCheck(data)
            data = self.SendAction(data)
            self.client.send(data)
            self.client.close()
        self.s.close()

class StatusServer(BasicServerSocket):
    STATUS = 'READY'
    def __init__(self, port, block_size, myname=None):
        self.port = port
        self.block_size = block_size
        BasicServerSocket.__init__(self, port, block_size, myname)
        self.setDaemon(True)

    def SendAction(self, data):
        if data == 'P_RESET_STATUS':
            self.STATUS = 'READY'
        data = self.STATUS
        return data

class BasicServer(BasicServerSocket):
    PROTOCOL = None
    debug = True
    RESULT = None
    BASIC_COMMAND_LIST = ('P_GETDATA', 'P_STORE_DATA', 'P_NONE')
    COMMAND_LIST = ()
    status_server = None
    STATUS = 'READY'

    def __init__(self, port, block_size, status_server=None, myname=None):
        self.port = port
        self.block_size = block_size
        ##  self.COMMAND_LIST = ()
        if status_server != None:
            self.status_server = status_server
        self.PROTOCOL = {}
        self.BuildProtocolTable()
        BasicServerSocket.__init__(self, port, block_size, myname)

    def setStatus(self,status):
        self.STATUS = status
        if self.status_server != None:
            self.status_server.STATUS = status

    def SendAction(self, data):
        data = data.split(',')
        if data[0] in list(self.PROTOCOL.keys()):
            print(self.client_address[0] + ', ' + time.strftime('%H:%M:%S') +\
            ', ' + self.getName() + ', EXECUTE ' + str(data).replace(',', ''))
            try:
                if len(data) > 1:
                    data = str(self.PROTOCOL[data[0]](data[1:]))
                else:
                    data = str(self.PROTOCOL[data[0]]())
            except Exception as ex:
                print('ProcessException', ex)
                data = 'False'
        else:
            print(self.client_address[0] + ', ' + time.strftime('%H:%M:%S') +\
            ', ' + self.getName() + ', UNKNOWN ' + str(data).replace(',', ''))
            data = 'False'
        return data

    def P_GETDATA(self, *args):
        self.setStatus('SENDING_DATA')
        F = io.StringIO()
        pickle.dump(self.RESULT, F, PICKLE_PROTOCOL)
        F.seek(0)
        data = 'OK'
        while data != '':
            data = F.read(self.block_size)
            self.client.send(data)
        self.setStatus('READY')
        print(octopussy)
        return True

    def P_STORE_DATA(self, *args):
        global HOSTNAME
        G = open(HOSTNAME + '_data.bin','wb')
        pickle.dump(self.RESULT, G, PICKLE_PROTOCOL)
        G.flush()
        G.close()
        return True

    def P_NONE(self, *args):
        return True

    def BuildProtocolTable(self):
        for cmd in self.BASIC_COMMAND_LIST:
            self.PROTOCOL.setdefault(cmd,getattr(self,cmd))
        for cmd in self.COMMAND_LIST:
            self.PROTOCOL.setdefault(cmd,getattr(self,cmd))

class ModelFileServer(BasicServerSocket):
    model_file = None
    model_file_name = 'None'
    model_directory = None

    def __init__(self, port, block_size, myname=None):
        self.port = port
        self.block_size = block_size
        BasicServerSocket.__init__(self, port, block_size, myname)
        self.setDaemon(True)

    def ReadFile(self, model_file_name, model_directory=None):
        if self.model_directory != None and model_directory == None:
            model_directory = self.model_directory
        fullP = os.path.join(model_directory,model_file_name)
        if os.path.exists(fullP):
            self.model_file = open(fullP,'r')
            self.model_directory = model_directory
            self.model_file_name = model_file_name
            return True
        else:
            return False

    def ListenLoop(self):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.bind(('', self.port))
        self.s.listen(self.backlog)

        while self.server_active:
            self.client, self.client_address = self.s.accept()
            data = self.client.recv(self.block_size)
            data = self.KillCheck(data)
            if self.RequestLogOn: self.RequestLog(data)
            if data == 'GET':
                self.model_file.seek(0)
                while data != '':
                    data = self.model_file.read(self.block_size)
                    self.client.send(data)
            elif data == 'LIST':
                F = io.StringIO()
                data = os.listdir(self.model_directory)
                pickle.dump(data,F,PICKLE_PROTOCOL)
                F.seek(0)
                while data != '':
                    data = F.read(self.block_size)
                    self.client.send(data)
            elif data[:4] == 'LOAD':
                data = data.split(',')[1]
                data = str(self.ReadFile(data))
            else:
                data = 'False'
                self.client.send(data)
            self.client.close()
        self.s.close()

class ServerStatusCheck(Thread):
    """
    Creates a Thread that polls [servers] on port every interval seconds
    for their current status. This is collected in a current_status list
    as (server,status) tuples
    """
    servers = None
    port = None
    block_size = None
    interval = None
    current_status = None
    go = True

    def __init__(self, servers, port, block_size, interval=120, myname=None):
        Thread.__init__(self)
        self.servers = servers
        self.port = port
        self.block_size = block_size
        self.interval = interval
        self.current_status = []
        self.setDaemon(True)

    def run(self):
        while self.go:
            self.current_status = []
            for s in self.servers:
                self.current_status.append(self.PollServer(s))
            self.PrintStatus()
            time.sleep(self.interval)

    def PollServer(self, server):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.connect((server, self.port))
        self.s.send('STATUS')
        data = self.s.recv(self.block_size)
        self.s.close()
        return (server, data)

    def PrintStatus(self):
        print('\n********** Status report %s************\n*' % time.strftime('%H:%M:%S'))
        for s in self.current_status:
            print('*', s[0], s[1])
        print('*\n*********************************************\n')

class TentacleScanner:
    def __init__(self,servers):
        self.servers = servers
        self.servers_ready = []
        self.servers_busy = []
        self.servers_dead = []
        self.feedback = []
        self.feedback_history = []

    def scan(self):
        self.feedback = []
        self.servers_ready = []
        self.servers_busy = []
        self.servers_dead = []

        for server in self.servers:
            try:
                print('Tentacle scanner is trying server:', server)
                client = SimpleClient(server, STATUS_PORT, BLOCK_SIZE)
                client.timeout = 5
                client.send(['STATUS'])
                print('Response:', client.response)
                self.feedback.append((client.server,client.response[0]))
            except Exception as ex:
                print(ex)
                self.feedback.append((client.server,'FAILED'))
        self.feedback_history.append(self.feedback)

        ##  print self.feedback, '\n'

        for sv in self.feedback:
            if sv[1] == 'FAILED':
                self.servers_dead.append(sv[0])
            elif sv[1] == 'READY':
                self.servers_ready.append(sv[0])
            else:
                self.servers_busy.append(sv[0])
        print('\nready:\n%s \n' % self.servers_ready)
        print('busy:\n%s \n'  % self.servers_busy)
        print('dead:\n%s \n'  % self.servers_dead)

    def getAvailableServers(self):
        self.scan()
        return self.servers_ready

    def getWorkingServers(self):
        self.scan()
        return self.servers_busy

    def getActiveServers(self):
        self.scan()
        return self.servers_ready + self.servers_busy

    def getDeadServers(self):
        self.scan()
        return self.servers_dead

class ServerListLoader:
    file_name = 'server_list'
    directory_name = os.path.dirname(os.path.abspath(os.sys.argv[0]))
    server_list = None

    def __init__(self):
        self.server_list = []

    def ReadFile(self, file_name=None, directory_name=None):
        if file_name != None:
            self.file_name = file_name
        if directory_name != None:
            self.directory_name = directory_name

        self.server_list = []
        try:
            sFile = open(os.path.join(self.directory_name, self.file_name),'r')
            for l in sFile:
                l = l.strip()
                l = l.strip('\n')
                l = l.strip('\r')
                l = l.strip('\r\n')
                if l == '':
                    pass
                elif l[0] == '#':
                    pass
                else:
                    self.server_list.append(l)
            sFile.close()
            return self.server_list
        except Exception as ex:
            print(ex)
            print('Cannot find \'server_list\' file in current directory: %s' % self.directory_name)
            print('This is a fatal error please create this file with server names, one per line')
            return []
