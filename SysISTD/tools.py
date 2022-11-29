'''
Created on 24 Apr 2020

@author: Francois STÜDER
'''
import os
import shutil
from threading import Thread
from time import sleep

def test_file_can_open(file,errorMessage):
    try:
        with open(file):
            pass
    except IOError:
        raise Exception(errorMessage)

def make_dirs(d):
    try:
        os.makedirs(d)
        return True
    except:
        return False
    
def remove_dir(dir_path,printError=True):
    try:
        shutil.rmtree(dir_path)
    except OSError as e:
        if printError:
            print("Error: %s : %s" % (dir_path, e.strerror))
        
def move_dir(Dir,Folder):
    try:
        shutil.move(Dir,Folder)
        return True,None
    except Exception as e:
        return False,e

def init_list_thread(ParrallelTask):
    '''Initialisation Threads'''
    ParrallelTaskList = []
    for t in range(ParrallelTask):
        ParrallelTaskList.append(initiate())
        ParrallelTaskList[-1].start()
    return ParrallelTaskList

def auto_add_thread(ListThread,WaitTime,job):
    while True:
        t = 0
        for busyThread in ListThread:
            if busyThread.is_alive() == False:
                ListThread[t] = job
                job.start()
                return True
            else:
                t = t+1
    sleep(WaitTime)

def wait_end_thread_list(ParrallelTaskList):
    for OneThread in ParrallelTaskList:
        if OneThread.is_alive() == True:
            OneThread.join()

class initiate(Thread):

    def __init__(self):
        Thread.__init__(self)

    def run(self):
        return True

core = 1
tmp = "/tmp"