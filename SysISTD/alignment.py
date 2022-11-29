'''
Created on 18 Nov 2022

@author: Francois STÜDER
'''

from io import TextIOWrapper
#from multiprocessing import Pool
from os import path,remove
from subprocess import Popen, PIPE
#from threading import Thread

from SysISTD import tools

cmdAlignment = "bowtie2 --reorder --very-sensitive -p __core__ -x __genome__ -1 __file1__ -2 __file2__ -S __outputSam__"
samDir = "./"

def alignment_generate_sam_coordinate(file1, file2, genomePath):
    """
    Alignment : flexible cmd with config file
    Keywords :
    __file1__ : fastq file 1
    __file2__ : fastq file 2
    __core__ : Number core use
    __genome__ : Genome path
    """
    cmd_alignment = []
    global samDir
    samFileOut = path.join(samDir, "all.sam")
    global cmdAlignment
    for partCMD in cmdAlignment.split(' '):
        if partCMD == '__file1__':
            cmd_alignment.append(file1)
        elif partCMD == '__file2__':
            cmd_alignment.append(file2)
        elif partCMD == '__core__':
            cmd_alignment.append(str(tools.core))
        elif partCMD == '__genome__':
            cmd_alignment.append(genomePath)
        elif partCMD == '__outputSam__':
            cmd_alignment.append(samFileOut)
        else:
            cmd_alignment.append(partCMD)
    alignRun = Popen(cmd_alignment, stdout=PIPE, stderr=PIPE)
    readMapped = 0
    header = ""
    coordinateAlreadyWrited = set()
    stdout, stderror =  alignRun.communicate()
    
    with open(path.join(samDir, f'alignment.err'), 'w') as logerr:
        logerr.write(stderror.decode("utf-8"))
    
    samFileOutHandler = open(samFileOut, 'r')
    line = next(samFileOutHandler)
    while line.startswith('@'):
        header += line
        line = next(samFileOutHandler)
    
    cutedLine = line.split('\t')
    previousBarcode = cutedLine[0].split("/")[0].split(":")[-1]
    samOut = open(path.join(samDir, f'{previousBarcode}.sam'), 'a')
    samOut.write(header)
    coordinateAlreadyWrited.add(previousBarcode)
    for line in samFileOutHandler:
        cutedLine = line.split('\t')
        barcode = cutedLine[0].split("/")[0].split(":")[-1]
        if barcode != previousBarcode:
            previousBarcode = barcode
            samOut.close()
            samOut = open(path.join(samDir, f'{barcode}.sam'), 'a')
            if barcode not in coordinateAlreadyWrited:
                samOut.write(header)
                coordinateAlreadyWrited.add(barcode)
        samOut.write(line)
        if cutedLine[2] != '*':
            readMapped += 0.5
    samFileOutHandler.close()
    
    remove(samFileOut)
    return readMapped

'''def alignment_generate_sam_coordinate(file1, file2, genomePath):
    """
    Alignment : flexible cmd with config file
    Keywords :
    __file1__ : fastq file 1
    __file2__ : fastq file 2
    __core__ : Number core use
    __genome__ : Genome path
    """
    cmd_alignment = []
    global cmdAlignment
    for partCMD in cmdAlignment.split(' '):
        if partCMD == '__file1__':
            cmd_alignment.append(file1)
        elif partCMD == '__file2__':
            cmd_alignment.append(file2)
        elif partCMD == '__core__':
            cmd_alignment.append(str(tools.core))
        elif partCMD == '__genome__':
            cmd_alignment.append(genomePath)
        else:
            cmd_alignment.append(partCMD)
    alignRun = Popen(cmd_alignment, stdout=PIPE, stderr=PIPE)

    readMapped = 0
    header = ""
    coordinateAlreadyWrited = set()
    global samDir
    writePool = Pool(processes=1)
    stdoutFileHandler = TextIOWrapper(alignRun.stdout, encoding="utf-8")
    line = next(stdoutFileHandler)
    while line.startswith('@'):
        header += line
        line = next(stdoutFileHandler)
    oldreadLines = get_lines(stdoutFileHandler)
    oldreadLines.start()
    writePool.apply_async(write_sam, ([line], header, coordinateAlreadyWrited, readMapped))
    readLines = get_lines(stdoutFileHandler)
    oldreadLines.join()
    readLines.start()
    while oldreadLines.alreadyReadable:
        writePool.apply_async(write_sam, (oldreadLines.groupLine, header, coordinateAlreadyWrited, readMapped))
        oldreadLines = readLines
        readLines = get_lines(stdoutFileHandler)
        oldreadLines.join()
        readLines.start()
    lastWrite = writePool.apply_async(write_sam, (oldreadLines.groupLine, header, coordinateAlreadyWrited, readMapped))
    with open(path.join(samDir, f'alignment.err'), 'w') as logerr:
        for line in TextIOWrapper(alignRun.stderr, encoding="utf-8"):
            logerr.write(line)
    lastWrite.get()
    #if alignRun.returncode != 0:
    #    raise Exception(f'Error with alignment !')

    return readMapped

class get_lines(Thread):
    """
    classdocs
    """


    def __init__(self, fileHandler):
        """
        Constructor
        """
        Thread.__init__(self)
        self.fileHandler = fileHandler
        self.groupLine = []
        self.alreadyReadable = True

    def run(self):
        try:
            for i in range(10000):
                line = next(self.fileHandler)
                self.groupLine.append(line)
        except:
            self.alreadyReadable = False

def write_sam(groupLine, header, coordinateAlreadyWrited, readMapped):
    if len(groupLine) > 0:
        cutedLine = groupLine.pop(0).split('\t')
        previousBarcode = cutedLine[0].split("/")[0].split(":")[-1]
        samOut = open(path.join(samDir, f'{previousBarcode}.sam'), 'a')
        if previousBarcode not in coordinateAlreadyWrited:
            samOut.write(header)
            coordinateAlreadyWrited.add(previousBarcode)
        for line in groupLine:
            cutedLine = line.split('\t')
            barcode = cutedLine[0].split("/")[0].split(":")[-1]
            if barcode != previousBarcode:
                previousBarcode = barcode
                samOut.close()
                samOut = open(path.join(samDir, f'{barcode}.sam'), 'a')
                if barcode not in coordinateAlreadyWrited:
                    samOut.write(header)
                    coordinateAlreadyWrited.add(barcode)
            samOut.write(line)
            if cutedLine[2] != '*':
                readMapped += 0.5
        
        samOut.close()'''