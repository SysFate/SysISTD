'''
Created on 17 Nov 2022

@author: Francois STÜDER
'''

import regex

from os import path
from threading import Thread,RLock
from time import time

guibsonSeq = {'+':'XXX','-':'XXX'}
polyT = {"+": "XXX", "-": "XXX", "s+": "XXX", "s-": "XXX"}
polyTPair = {'+':'XXX','-':'XXX'}
barcodesX = {}
barcodesY = {}
lengthBC = 0
calculPool = None
polyTSearch = True
readPas = 3000
UMISearch = True
trimSequence = 0
lengthUMI = 9
fastqDir = "./"
outBadSeqFolder = "./"
goodRead = RLock()
verrouStatsOut = RLock()
verrouInsertUMIFind = RLock()
writeGoodRead = RLock()
verrouInsertInQueueWriteGoodRead = RLock()
queueWriteGoodRead = []

class AnnalyseRecord(Thread):
    '''
    classdocs
    '''


    def __init__(self, listRecord,dicCountFindReadByCoordinate,statsOut,UMIStats,statTimeTake):
        '''
        Constructor
        '''
        Thread.__init__(self)
        self.listRecord = listRecord
        self.dicCountFindReadByCoordinate = dicCountFindReadByCoordinate
        self.statsOut = statsOut
        self.UMIStats = UMIStats
        self.statTimeTake = statTimeTake

    def run(self):
        beginWait = time()
        global guibsonSeq
        global polyT
        global polyTPair
        global barcodesX
        global barcodesY
        global lengthBC
        global calculPool
        global polyTSearch
        global readPas
        global UMISearch
        listJob = []
        i = 0
        while len(self.listRecord) > i:
            iEnd = i+readPas
            if polyTSearch:
                job = calculPool.apply_async(extracts_informations_reads, (self.listRecord[i:iEnd],guibsonSeq,polyT,polyTPair,barcodesX,barcodesY,lengthBC))
            else:
                job = calculPool.apply_async(extracts_informations_reads_without_polyT, (self.listRecord[i:iEnd],guibsonSeq,barcodesX,barcodesY,lengthBC))
            listJob.append(job)
            i = iEnd
        
        self.listRecord = None
        localDicFastqCoordinate = {}
        localDicCountFindReadByCoordinate = {}
        localUMIStats = {}
        statsOutLocal = {-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,'total_read_count':0,'goodPolyT':0,'goodG':0,'reads_good_bacode':0}
        for job in listJob:
            results = job.get()
            for result in results:
                statsOutLocal['total_read_count'] += 1
                statsOutLocal[result[0]] = statsOutLocal[result[0]] + 1
                if result[0] == 1:
                    statsOutLocal['reads_good_bacode'] += 1
                    if result[3] not in localDicCountFindReadByCoordinate:
                        localDicCountFindReadByCoordinate[result[3]] = 1
                    else:
                        localDicCountFindReadByCoordinate[result[3]] += 1
                    key1 = f'{result[3]}_1'
                    key2 = f'{result[3]}_2'
                    if key1 not in localDicFastqCoordinate:
                        localDicFastqCoordinate[key1] = result[1].format('fastq')
                    else:
                        localDicFastqCoordinate[key1] = "{}{}".format(localDicFastqCoordinate[key1], result[1].format('fastq'))
                    if key2 not in localDicFastqCoordinate:
                        localDicFastqCoordinate[key2] = result[2].format('fastq')
                    else:
                        localDicFastqCoordinate[key2] = "{}{}".format(localDicFastqCoordinate[key2], result[2].format('fastq'))
                    if UMISearch:
                        if result[3] not in localUMIStats:
                            localUMIStats[result[3]] = {}
                        for barcodePos,seqsUMI in result[4].items():
                            if barcodePos not in localUMIStats[result[3]]:
                                localUMIStats[result[3]][barcodePos] = {}
                            for seqUMI in seqsUMI:
                                if seqUMI not in localUMIStats[result[3]][barcodePos]:
                                    localUMIStats[result[3]][barcodePos][seqUMI] = 1
                                else:
                                    localUMIStats[result[3]][barcodePos][seqUMI] = localUMIStats[result[3]][barcodePos][seqUMI] + 1
                elif result[0] == -4:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == -3:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == -2:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == -1:
                    statsOutLocal['goodG'] += 1
                elif result[0] == 0:
                    pass
                elif result[0] == 2:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 3:
                    statsOutLocal['goodG'] += 1
                elif result[0] == 4:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 5:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 6:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 7:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 8:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 9:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 10:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 11:
                    statsOutLocal['goodPolyT'] += 1
                elif result[0] == 12:
                    statsOutLocal['goodPolyT'] += 1
        
        listJob = None
        beginWaitLock = time()
        with goodRead:
            self.statTimeTake['get_listRead_lock'] += (time()-beginWaitLock)
            for key,value in localDicCountFindReadByCoordinate.items():
                if key not in self.dicCountFindReadByCoordinate:
                    self.dicCountFindReadByCoordinate[key] = value
                else:
                    self.dicCountFindReadByCoordinate[key] += value
        localDicCountFindReadByCoordinate = None

        beginWaitLock = time()
        with verrouStatsOut:
            self.statTimeTake['get_stats_lock'] += (time()-beginWaitLock)
            for key,value in statsOutLocal.items():
                self.statsOut[key] = self.statsOut[key] + value
        statsOutLocal = None

        beginWaitLock = time()
        with verrouInsertUMIFind:
            self.statTimeTake['get_UMIStats_lock'] += (time()-beginWaitLock)
            for key,value in localUMIStats.items():
                if key not in self.UMIStats:
                    self.UMIStats[key] = value
                else:
                    for barcodePos,seqsUMI in value.items():
                        if barcodePos not in self.UMIStats[key]:
                            self.UMIStats[key][barcodePos] = seqsUMI
                        else:
                            for seqUMI,counts in seqsUMI.items():
                                if seqUMI not in self.UMIStats[key][barcodePos]:
                                    self.UMIStats[key][barcodePos][seqUMI] = counts
                                else:
                                    self.UMIStats[key][barcodePos][seqUMI] = self.UMIStats[key][barcodePos][seqUMI] + counts

        global queueWriteGoodRead
        beginWaitLock = time()
        if writeGoodRead.acquire(False):
            verrouInsertInQueueWriteGoodRead.acquire()
            self.statTimeTake['get_queue_lock'] += (time()-beginWaitLock)
            queueWriteGoodReadLocal = queueWriteGoodRead
            queueWriteGoodRead = []
            verrouInsertInQueueWriteGoodRead.release()
            merge_dic_fastq_coordinate(localDicFastqCoordinate, queueWriteGoodReadLocal)
            write_coordinate_fastq(localDicFastqCoordinate)
            writeGoodRead.release()
            queueWriteGoodReadLocal = None
        else:
            verrouInsertInQueueWriteGoodRead.acquire()
            self.statTimeTake['get_queue_lock'] += (time()-beginWaitLock)
            if len(queueWriteGoodRead) >= 10:
                queueWriteGoodReadLocal = queueWriteGoodRead
                queueWriteGoodRead = []
                verrouInsertInQueueWriteGoodRead.release()
                merge_dic_fastq_coordinate(localDicFastqCoordinate, queueWriteGoodReadLocal)
                queueWriteGoodReadLocal = None
                verrouInsertInQueueWriteGoodRead.acquire()

            queueWriteGoodRead.append(localDicFastqCoordinate)
            verrouInsertInQueueWriteGoodRead.release()

        localUMIStats = None
        localDicFastqCoordinate = None
        value = None
        self.statTimeTake['all_analyse'] += (time()-beginWait)

""" This function will insert all others dictionar into the first give in argument.
It is design to merge coordinate fastq to create one
"""
def merge_dic_fastq_coordinate(mergedDicFastqCoordinate,listOfOtherDic):
    for otherDic in listOfOtherDic:
        for key, value in otherDic.items():
            if key in mergedDicFastqCoordinate:
                mergedDicFastqCoordinate[key] = "".join([mergedDicFastqCoordinate[key], value])
            else:
                mergedDicFastqCoordinate[key] = value

def write_coordinate_fastq(dicFastqCoordinateToWrite):
    global fastqDir
    for key, value in dicFastqCoordinateToWrite.items():
        fastqCoordinate = open(path.join(fastqDir, f'{key}.fastq'), 'a')
        fastqCoordinate.write(value)
        fastqCoordinate.close()

def extracts_informations_reads(listRecords,guibsonSeq,polyT,polyTPair,barcodesX,barcodesY,lengthBC):
    dataExtracts = []
    for records in listRecords:
        dataExtracts.append(extract_good_reads(records[0],records[1],guibsonSeq,polyT,polyTPair,barcodesX,barcodesY,lengthBC))
    listRecords = None
    return dataExtracts
        
def extract_good_reads(record1,record2,guibsonSeq,polyT,polyTPair,barcodesX,barcodesY,lengthBC):
    """
    Permet de determier si notre sequence est conforme aux attentes.
        -4 : Two good construct in the same positions
        -3 : Too many good construct (more than 2)
        -2 : Too many good construct and bad construct
        -1 : Only bad construct
        0 : No Guibson
        1 : Good
        2 : Barcode strand problem
        3 : No polyT found in the right position
        4 : No Barcode found
        5 : Multiple barcode found
    This list are when we have two good contruct
        6 : Barcorde not the same
        7 : insertitude in multiple barcode on one read after validation
        8 : min one read have a recode error 2 and the rest can't validate a group of barcode
        9 : The rest
    For the test of barcode:
        10 : Multiple match BCX
        11 : Multiple match BCY
        12 : Multiple match BCX and BCY
    Extrait aussi les barcode ainsi que les bonnes sequences
    """
    # On recherche Guibson. Ainsi seul les reads contenant guibson sont accepte.
    find1P = extact_positions(guibsonSeq["+"], str(record1.seq))
    find1N = extact_positions(guibsonSeq["-"], str(record1.seq))
    find2P = extact_positions(guibsonSeq["+"], str(record2.seq))
    find2N = extact_positions(guibsonSeq["-"], str(record2.seq))
    
    findPolyT1P = None
    findPolyT2P = None
    findPolyT1N = None
    findPolyT2N = None
    if find1P != None:
        findPolyT1P = extact_positions(polyT["+"], str(record1.seq))
    if find2P != None:
        findPolyT2P = extact_positions(polyT["+"], str(record2.seq))
    if find1N != None:
        findPolyT1N = extact_positions(polyT["-"], str(record1.seq))
    if find2N != None:
        findPolyT2N = extact_positions(polyT["-"], str(record2.seq))
    
    # On va maintenant tester si il n'y a que 1 combinaison qui fonctionne sinon cela est un artefact
    goodComb = (findPolyT1P != None) + (findPolyT2P != None) + (findPolyT1N != None) + (findPolyT2N != None)
    
    if goodComb == 1:
        if findPolyT1P != None:
            findPolyT = findPolyT1P
            if (find1P[0]-lengthBC) < 0:
                XBC,numberXBC = extract_barcode(barcodesX,"0"*(lengthBC-find1P[0])+str(record1[0:find1P[0]].seq))
            else:
                XBC,numberXBC = extract_barcode(barcodesX,str(record1[find1P[0]-lengthBC:find1P[0]].seq))
            YBC,numberYBC = extract_barcode(barcodesY,str(record1[find1P[1]:find1P[1]+lengthBC].seq)+"0"*((find1P[1]+lengthBC)-len(record1)))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"+") #record1
            if readCode == 1:
                recordCut1 = record1[findPolyT[1]:]
                if len(recordCut1) == 0:
                    if len(recordCut1) == 0:
                        recordCut1 = record1[-1:] # Eviter les erreurs avec des sequences vide
                    findPolyTPair = extact_positions(polyTPair["-"], str(record2.seq))
                    if findPolyTPair != None:
                        recordCut2 = record2[:findPolyTPair[0]]
                    else:
                        recordCut2 = record2
                else:
                    recordCut2 = record2
                UMI = extract_UMI(record1,find1P,"+",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut1, recordCut2 = trim_sequences_two_records(
                    recordCut1, recordCut2, 'findPolyTPair' in locals() and findPolyTPair != None)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif findPolyT2P != None:
            findPolyT = findPolyT2P
            if (find2P[0]-lengthBC) < 0:
                XBC,numberXBC = extract_barcode(barcodesX,"0"*(lengthBC-find2P[0])+str(record2[0:find2P[0]].seq))
            else:
                XBC,numberXBC = extract_barcode(barcodesX,str(record2[find2P[0]-lengthBC:find2P[0]].seq))
            YBC,numberYBC = extract_barcode(barcodesY,str(record2[find2P[1]:find2P[1]+lengthBC].seq)+"0"*((find2P[1]+lengthBC)-len(record2)))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"+") #record2
            if readCode == 1:
                recordCut2 = record2[findPolyT[1]:]
                if len(recordCut2) == 0:
                    recordCut2 = record2[-1:] # Eviter les erreurs avec des sequences vide
                    findPolyTPair = extact_positions(polyTPair["-"], str(record1.seq))
                    if findPolyTPair != None:
                        recordCut1 = record1[:findPolyTPair[0]]
                    else:
                        recordCut1 = record1
                else:
                    recordCut1 = record1
                UMI = extract_UMI(record2,find2P,"+",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut2, recordCut1 = trim_sequences_two_records(
                    recordCut2, recordCut1, 'findPolyTPair' in locals() and findPolyTPair != None)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif findPolyT1N != None:
            findPolyT = findPolyT1N
            XBC,numberXBC = extract_barcode(barcodesX,str(record1[find1N[1]:find1N[1]+lengthBC].seq)+"0"*((find1N[1]+lengthBC)-len(record1)))
            if (find1N[0]-lengthBC) < 0:
                YBC,numberYBC = extract_barcode(barcodesY,"0"*(lengthBC-find1N[0])+str(record1[0:find1N[0]].seq))
            else:
                YBC,numberYBC = extract_barcode(barcodesY,str(record1[find1N[0]-lengthBC:find1N[0]].seq))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"-") #record1
            if readCode == 1:
                recordCut1 = record1[:findPolyT[0]]
                if len(recordCut1) == 0:
                    recordCut1 = record1[:1] # Eviter les erreurs avec des sequences vide
                    findPolyTPair = extact_positions(polyTPair["+"], str(record2.seq))
                    if findPolyTPair != None:
                        recordCut2 = record2[:findPolyTPair[0]]
                    else:
                        recordCut2 = record2
                else:
                    recordCut2 = record2
                UMI = extract_UMI(record1,find1N,"-",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut2 = trim_sequences_one_records(
                    recordCut2, 'findPolyTPair' in locals() and findPolyTPair != None)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif findPolyT2N != None:
            findPolyT = findPolyT2N
            XBC,numberXBC = extract_barcode(barcodesX,str(record2[find2N[1]:find2N[1]+lengthBC].seq)+"0"*((find2N[1]+lengthBC)-len(record2)))
            if (find2N[0]-lengthBC) < 0:
                YBC,numberYBC = extract_barcode(barcodesY,"0"*(lengthBC-find2N[0])+str(record2[0:find2N[0]].seq))
            else:
                YBC,numberYBC = extract_barcode(barcodesY,str(record2[find2N[0]-lengthBC:find2N[0]].seq))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"-") #record2
            if readCode == 1:
                recordCut2 = record2[:findPolyT[0]]
                if len(recordCut2) == 0:
                    recordCut2 = record2[:1] # Eviter les erreurs avec des sequences vide
                    findPolyTPair = extact_positions(polyTPair["+"], str(record1.seq))
                    if findPolyTPair != None:
                        recordCut1 = record1[:findPolyTPair[0]]
                    else:
                        recordCut1 = record1
                else:
                    recordCut1 = record1
                UMI = extract_UMI(record2,find2N,"-",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut1 = trim_sequences_one_records(
                    recordCut1, 'findPolyTPair' in locals() and findPolyTPair != None)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
    elif goodComb == 2:
        if findPolyT1P != None and findPolyT2N != None:
            readCode,recordCut1,recordCut2,coordinateName,UMI = search_P_N_record(find1P,find2N,findPolyT1P,findPolyT2N,record1,record2,barcodesX,barcodesY,lengthBC)
        elif findPolyT2P != None and findPolyT1N != None:
            readCode,recordCut2,recordCut1,coordinateName,UMI = search_P_N_record(find2P,find1N,findPolyT2P,findPolyT1N,record2,record1,barcodesX,barcodesY,lengthBC)
        else:
            readCode = -4
        if readCode == 1:
            return readCode,recordCut1,recordCut2,coordinateName,UMI
    else:
        guibFind = (find1N != None) + (find2N != None) + (find2P != None) + (find1P != None)
        if guibFind == 0:
            readCode = 0
        else:
            findPolyT1PS = None
            findPolyT2PS = None
            findPolyT1NS = None
            findPolyT2NS = None
            if find1P != None:
                findPolyT1PS = extact_positions(polyT["s+"], str(record1.seq))
            if find2P != None:
                findPolyT2PS = extact_positions(polyT["s+"], str(record2.seq))
            if find1N != None:
                findPolyT1NS = extact_positions(polyT["s-"], str(record1.seq))
            if find2N != None:
                findPolyT2NS = extact_positions(polyT["s-"], str(record2.seq))
            badComb = (findPolyT1PS != None) + (findPolyT2PS != None) + (findPolyT1NS != None) + (findPolyT2NS != None)
            if goodComb == 0:
                if badComb > 0:
                    readCode = -1
                else:
                    readCode = 3
            else:
                if badComb > 0:
                    readCode = -2
                else:
                    readCode = -3
    return readCode,record1,record2

def extracts_informations_reads_without_polyT(listRecords, guibsonSeq, barcodesX, barcodesY, lengthBC):
    dataExtracts = []
    for records in listRecords:
        dataExtracts.append(extract_good_reads_without_polyT(records[0],records[1],guibsonSeq,barcodesX,barcodesY,lengthBC))
    listRecords = None
    return dataExtracts

def extract_good_reads_without_polyT(record1,record2,guibsonSeq,barcodesX,barcodesY,lengthBC):
    """
    Permet de determier si notre sequence est conforme aux attentes.
        -4 : Two good construct in the same positions
        -3 : Too many good construct (more than 2)
        -2 : Too many good construct and bad construct
        -1 : Only bad construct
        0 : No Guibson
        1 : Good
        2 : Barcode strand problem
        3 : No polyT found in the right position
        4 : No Barcode found
        5 : Multiple barcode found
    This list are when we have two good contruct
        6 : Barcorde not the same
        7 : insertitude in multiple barcode on one read after validation
        8 : min one read have a recode error 2 and the rest can't validate a group of barcode
        9 : The rest
    For the test of barcode:
        10 : Multiple match BCX
        11 : Multiple match BCY
        12 : Multiple match BCX and BCY
    Extrait aussi les barcode ainsi que les bonnes sequences
    """
    # On recherche Guibson. Ainsi seul les reads contenant guibson sont accepte.
    find1P = extact_positions(guibsonSeq["+"], str(record1.seq))
    find1N = extact_positions(guibsonSeq["-"], str(record1.seq))
    find2P = extact_positions(guibsonSeq["+"], str(record2.seq))
    find2N = extact_positions(guibsonSeq["-"], str(record2.seq))
    
    goodComb = (find1N != None) + (find2N != None) + (find2P != None) + (find1P != None)
        
    if goodComb == 1:
        global lengthUMI
        if find1P != None:
            if (find1P[0]-lengthBC) < 0:
                XBC,numberXBC = extract_barcode(barcodesX,"0"*(lengthBC-find1P[0])+str(record1[0:find1P[0]].seq))
            else:
                XBC,numberXBC = extract_barcode(barcodesX,str(record1[find1P[0]-lengthBC:find1P[0]].seq))
            YBC,numberYBC = extract_barcode(barcodesY,str(record1[find1P[1]:find1P[1]+lengthBC].seq)+"0"*((find1P[1]+lengthBC)-len(record1)))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"+") #record1
            if readCode == 1:
                recordCut1 = record1[find1P[1]+lengthBC+lengthUMI:]
                if len(recordCut1) == 0:
                    if len(recordCut1) == 0:
                        recordCut1 = record1[-1:] # Eviter les erreurs avec des sequences vide
                recordCut2 = record2
                UMI = extract_UMI(record1,find1P,"+",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut1, recordCut2 = trim_sequences_two_records(recordCut1, recordCut2, False)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif find2P != None:
            if (find2P[0]-lengthBC) < 0:
                XBC,numberXBC = extract_barcode(barcodesX,"0"*(lengthBC-find2P[0])+str(record2[0:find2P[0]].seq))
            else:
                XBC,numberXBC = extract_barcode(barcodesX,str(record2[find2P[0]-lengthBC:find2P[0]].seq))
            YBC,numberYBC = extract_barcode(barcodesY,str(record2[find2P[1]:find2P[1]+lengthBC].seq)+"0"*((find2P[1]+lengthBC)-len(record2)))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"+") #record2
            if readCode == 1:
                recordCut2 = record2[find2P[1]+lengthBC+lengthUMI:]
                if len(recordCut2) == 0:
                    recordCut2 = record2[-1:] # Eviter les erreurs avec des sequences vide
                recordCut1 = record1
                UMI = extract_UMI(record2,find2P,"+",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut2, recordCut1 = trim_sequences_two_records(recordCut2, recordCut1, False)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif find1N != None:
            XBC,numberXBC = extract_barcode(barcodesX,str(record1[find1N[1]:find1N[1]+lengthBC].seq)+"0"*((find1N[1]+lengthBC)-len(record1)))
            if (find1N[0]-lengthBC) < 0:
                YBC,numberYBC = extract_barcode(barcodesY,"0"*(lengthBC-find1N[0])+str(record1[0:find1N[0]].seq))
            else:
                YBC,numberYBC = extract_barcode(barcodesY,str(record1[find1N[0]-lengthBC:find1N[0]].seq))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"-") #record1
            if readCode == 1:
                if (find1N[0]-lengthBC-lengthUMI) > 0:
                    recordCut1 = record1[:find1N[0]-lengthBC-lengthUMI]
                else:
                    recordCut1 = record1[:1] # Eviter les erreurs avec des sequences vide
                recordCut2 = record2
                UMI = extract_UMI(record1,find1N,"-",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut1, recordCut2 = trim_sequences_two_records(recordCut1, recordCut2, False)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
        elif find2N != None:
            XBC,numberXBC = extract_barcode(barcodesX,str(record2[find2N[1]:find2N[1]+lengthBC].seq)+"0"*((find2N[1]+lengthBC)-len(record2)))
            if (find2N[0]-lengthBC) < 0:
                YBC,numberYBC = extract_barcode(barcodesY,"0"*(lengthBC-find2N[0])+str(record2[0:find2N[0]].seq))
            else:
                YBC,numberYBC = extract_barcode(barcodesY,str(record2[find2N[0]-lengthBC:find2N[0]].seq))
            readCode,seqBC,coordinateName = extact_data_barcode(numberXBC,numberYBC,XBC,YBC,"-") #record2
            if readCode == 1:
                if (find2N[0]-lengthBC-lengthUMI) > 0:
                    recordCut2 = record2[:find2N[0]]
                else:
                    recordCut2 = record2[:1] # Eviter les erreurs avec des sequences vide
                recordCut1 = record1
                UMI = extract_UMI(record2,find2N,"-",lengthBC)
                rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI)
                recordCut2, recordCut1 = trim_sequences_two_records(recordCut2, recordCut1, False)
                return readCode,recordCut1,recordCut2,coordinateName,UMI
    elif goodComb == 2:
        if find1P != None and find2N != None:
            readCode,recordCut1,recordCut2,coordinateName,UMI = search_P_N_record_without_polyT(find1P,find2N,record1,record2,barcodesX,barcodesY,lengthBC)
        elif find2P != None and find1N != None:
            readCode,recordCut2,recordCut1,coordinateName,UMI = search_P_N_record_without_polyT(find2P,find1N,record2,record1,barcodesX,barcodesY,lengthBC)
        else:
            readCode = -4
        if readCode == 1:
            return readCode,recordCut1,recordCut2,coordinateName,UMI
    else:
        guibFind = (find1N != None) + (find2N != None) + (find2P != None) + (find1P != None)
        if guibFind == 0:
            readCode = 0
        else:
            readCode = 3
    return readCode,record1,record2

def extact_positions(querry,text):
    find = regex.search(querry,text,concurrent=True)
    if find != None and len(find.spans()) != 0:
        return find.spans()[0]
    else:
        None

def extract_barcode(barcodes,seq):
    # On recherche les barcodes dans la sequence donne
    # On revoie un dico de dico. La clef du premier est le barcode puis la clef puis suivant le strand enfin il y a les positions.
    results = {}
    numberBCResults = {}
    for barcode,querrys in barcodes.items():
        barcodeFind = False
        for strand,querry in querrys.items():
            find = regex.search(querry,seq,concurrent=True)
            if find != None and len(find.spans()) != 0:
                error = find.fuzzy_counts[0] + find.fuzzy_counts[1] + find.fuzzy_counts[2]
                if error not in results:
                    results[error] = {}
                    numberBCResults[error] = 0
                if barcode not in results[error]:
                    results[error][barcode] = {}
                    barcodeFind = True
                results[error][barcode][strand] = seq #find.spans()[0]
        if barcodeFind:
            numberBCResults[error] = numberBCResults[error] + 1
    if len(results) == 0:
        return None,0
    else:
        errorList = list(results.keys())
        return results[min(errorList)],numberBCResults[min(errorList)]
    
def extact_data_barcode(numberXBC,numberYBC,XBC,YBC,strand): #record
    if numberXBC == 1 and numberYBC == 1:
        for BCX,strandsX in XBC.items():
            for strandX, seqX in strandsX.items(): #coordinateX
                for BCY,strandsY in YBC.items():
                    for strandY, seqY in strandsY.items(): #coordinateY
                        # Ici seul les barcode ayant le sens du guibson sont accepte
                        if strandX == strand and strandY == strand:
                            coordinateName = f'{BCX}x{BCY}'
                            seqBC = {"x": seqX,"y": seqY} # str(record[coordinateX[0]:coordinateX[1]].seq) str(record[coordinateY[0]:coordinateY[1]].seq)
                            return 1,seqBC,coordinateName
        return 2,None,None
    elif numberXBC == 0 and numberYBC == 0:
        return 4,None,None
    elif numberXBC > 1 and numberYBC > 1:
        return 12,None,None
    elif numberXBC > 1:
        return 10,None,None
    elif numberYBC > 1:
        return 11,None,None
    else:
        return 5,None,None
    
def rename_record_barcoded(recordCut1,recordCut2,coordinateName,UMI):
    if 'x' in UMI:
        recordCut1.id = recordCut1.id[:-2]+":{}:{}:{}".format(UMI['x'][0],UMI['y'][0],coordinateName)+recordCut1.id[-2:]
        recordCut2.id = recordCut2.id[:-2]+":{}:{}:{}".format(UMI['x'][0],UMI['y'][0],coordinateName)+recordCut2.id[-2:]
    else:
        recordCut1.id = recordCut1.id[:-2]+":{}:{}:{}".format("",UMI['y'][0],coordinateName)+recordCut1.id[-2:]
        recordCut2.id = recordCut2.id[:-2]+":{}:{}:{}".format("",UMI['y'][0],coordinateName)+recordCut2.id[-2:]
    
def search_P_N_record(findP,findN,findPolyTP,findPolyTN,recordP,recordN,barcodesX,barcodesY,lengthBC):
    #Chercher les barcodes pour les deux
    if (findP[0]-lengthBC) < 0:
        XBC2,numberXBC2 = extract_barcode(barcodesX,"N"*(lengthBC-findP[0])+str(recordP[0:findP[0]].seq))
    else:
        XBC2,numberXBC2 = extract_barcode(barcodesX,str(recordP[findP[0]-lengthBC:findP[0]].seq))
    YBC2,numberYBC2 = extract_barcode(barcodesY,str(recordP[findP[1]:findP[1]+lengthBC].seq)+"N"*((findP[1]+lengthBC)-len(recordP)))
    
    XBC1,numberXBC1 = extract_barcode(barcodesX,str(recordN[findN[1]:findN[1]+lengthBC].seq)+"N"*((findN[1]+lengthBC)-len(recordN)))
    if (findN[0]-lengthBC) < 0 :
        YBC1,numberYBC1 = extract_barcode(barcodesY,"N"*(lengthBC-findN[0])+str(recordN[0:findN[0]].seq))
    else:
        YBC1,numberYBC1 = extract_barcode(barcodesY,str(recordN[findN[0]-lengthBC:findN[0]].seq))

    readCode2,seqBC2,coordinateName2 = extact_data_barcode(numberXBC2,numberYBC2,XBC2,YBC2,"+")
    readCode1,seqBC1,coordinateName1 = extact_data_barcode(numberXBC1,numberYBC1,XBC1,YBC1,"-")
    
    if readCode2 == 1 and readCode1 == 1:
        if coordinateName1 == coordinateName2:
            readCode = 1
            seqBC = seqBC1
            coordinateName = coordinateName1
            UMI = extract_UMI(recordP,findP,"+",lengthBC)
            for strand,listSeq in extract_UMI(recordN,findN,"-",lengthBC).items():
                if strand in UMI:
                    UMI[strand].extend(listSeq)
                else:
                    UMI[strand] = listSeq
        else:
            readCode = 6
    elif readCode2 == 1:
        readCode = 1
        seqBC = seqBC2
        coordinateName = coordinateName2
        UMI = extract_UMI(recordP,findP,"+",lengthBC)
    elif readCode1 == 1:
        readCode = 1
        seqBC = seqBC1
        coordinateName = coordinateName1
        UMI = extract_UMI(recordN,findN,"-",lengthBC)
    elif numberXBC1 > 1 or numberYBC1 > 1 or numberXBC2 > 1 or numberYBC2 > 1:
        readCode = 7
    elif readCode1 == 2 or readCode2 == 2:
        readCode = 8
    else:
        readCode = 9
    
    if readCode == 1:
        recordCutN = recordN[:findPolyTN[0]]
        if len(recordCutN) == 0:
            recordCutN = recordN[:1] # Eviter les erreurs avec des sequences vide    
        recordCutP = recordP[findPolyTP[1]:]
        if len(recordCutP) == 0:
            recordCutP = recordP[-1:] # Eviter les erreurs avec des sequences vide
        rename_record_barcoded(recordCutP,recordCutN,coordinateName,UMI)
        return readCode,recordCutP,recordCutN,coordinateName,UMI
    else:
        return readCode,None,None,None,None

def search_P_N_record_without_polyT(findP,findN,recordP,recordN,barcodesX,barcodesY,lengthBC):
    #Chercher les barcodes pour les deux
    if (findP[0]-lengthBC) < 0:
        XBC2,numberXBC2 = extract_barcode(barcodesX,"N"*(lengthBC-findP[0])+str(recordP[0:findP[0]].seq))
    else:
        XBC2,numberXBC2 = extract_barcode(barcodesX,str(recordP[findP[0]-lengthBC:findP[0]].seq))
    YBC2,numberYBC2 = extract_barcode(barcodesY,str(recordP[findP[1]:findP[1]+lengthBC].seq)+"N"*((findP[1]+lengthBC)-len(recordP)))
    
    XBC1,numberXBC1 = extract_barcode(barcodesX,str(recordN[findN[1]:findN[1]+lengthBC].seq)+"N"*((findN[1]+lengthBC)-len(recordN)))
    if (findN[0]-lengthBC) < 0 :
        YBC1,numberYBC1 = extract_barcode(barcodesY,"N"*(lengthBC-findN[0])+str(recordN[0:findN[0]].seq))
    else:
        YBC1,numberYBC1 = extract_barcode(barcodesY,str(recordN[findN[0]-lengthBC:findN[0]].seq))

    readCode2,seqBC2,coordinateName2 = extact_data_barcode(numberXBC2,numberYBC2,XBC2,YBC2,"+")
    readCode1,seqBC1,coordinateName1 = extact_data_barcode(numberXBC1,numberYBC1,XBC1,YBC1,"-")
    
    if readCode2 == 1 and readCode1 == 1:
        if coordinateName1 == coordinateName2:
            readCode = 1
            seqBC = seqBC1
            coordinateName = coordinateName1
            UMI = extract_UMI(recordP,findP,"+",lengthBC)
            for strand,listSeq in extract_UMI(recordN,findN,"-",lengthBC).items():
                if strand in UMI:
                    UMI[strand].extend(listSeq)
                else:
                    UMI[strand] = listSeq
        else:
            readCode = 6
    elif readCode2 == 1:
        readCode = 1
        seqBC = seqBC2
        coordinateName = coordinateName2
        UMI = extract_UMI(recordP,findP,"+",lengthBC)
    elif readCode1 == 1:
        readCode = 1
        seqBC = seqBC1
        coordinateName = coordinateName1
        UMI = extract_UMI(recordN,findN,"-",lengthBC)
    elif numberXBC1 > 1 or numberYBC1 > 1 or numberXBC2 > 1 or numberYBC2 > 1:
        readCode = 7
    elif readCode1 == 2 or readCode2 == 2:
        readCode = 8
    else:
        readCode = 9
    
    if readCode == 1:
        global lengthUMI
        if findN[0]-lengthBC-lengthUMI > 0:
            recordCutN = recordN[:findN[0]-lengthBC-lengthUMI]
        else:
            recordCutN = recordN[:1] # Eviter les erreurs avec des sequences vide    
        recordCutP = recordP[findP[1]+lengthBC+lengthUMI:]
        if len(recordCutP) == 0:
            recordCutP = recordP[-1:] # Eviter les erreurs avec des sequences vide
        rename_record_barcoded(recordCutP,recordCutN,coordinateName,UMI)
        return readCode,recordCutP,recordCutN,coordinateName,UMI
    else:
        return readCode,None,None,None,None

def extract_UMI(record,find,strand,lengthBC):
    global lengthUMI
    UMI = {}
    if strand == '+':
        if (find[0]-lengthBC-lengthUMI) >= 0:
            UMI['x'] = [str(record[(find[0]-lengthBC-lengthUMI):(find[0]-lengthBC)].seq)]
        UMI['y'] = [str(record[(find[1]+lengthBC):(find[1]+lengthBC+lengthUMI)].seq)]
    elif strand == '-':
        if (find[1]+lengthBC+lengthUMI) <= len(record):
            UMI['x'] = [str(record[(find[1]+lengthBC):(find[1]+lengthBC+lengthUMI)].seq.reverse_complement())]
        UMI['y'] = [str(record[(find[0]-lengthBC-lengthUMI):(find[0]-lengthBC)].seq.reverse_complement())]
        
    return UMI

def trim_sequences_two_records(recordCut1, recordCut2, findPolyTPairIsNone):
    global trimSequence
    if isinstance(trimSequence, int) and trimSequence > 0:
        if len(recordCut1) <= trimSequence:
            recordCut1 = recordCut1[-1:]
        else:
            recordCut1 = recordCut1[:len(recordCut1)-trimSequence]
        if findPolyTPairIsNone:
            pass
        elif len(recordCut2) <= trimSequence:
            recordCut2 = recordCut2[-1:]
        else:
            recordCut2 = recordCut2[:len(recordCut2)-trimSequence]
    elif isinstance(trimSequence, str):
        pass
    return recordCut1, recordCut2

def trim_sequences_one_records(recordCut, findPolyTPairIsNone):
    global trimSequence
    if isinstance(trimSequence, int) and trimSequence > 0:
        if findPolyTPairIsNone:
            pass
        elif len(recordCut) <= trimSequence:
            recordCut = recordCut[-1:]
        else:
            recordCut = recordCut[:len(recordCut)-trimSequence]
    elif isinstance(trimSequence, str):
        pass
    return recordCut

def get_stats_code_sequences():
    return """
    Permet de determier si notre sequence est conforme aux attentes.
        -4 : Two good construct in the same positions
        -3 : Too many good construct (more than 2)
        -2 : Too many good construct and bad construct
        -1 : Only bad construct
        0 : No Guibson
        1 : Good
        2 : Barcode strand problem
        3 : No polyT found in the right position
        4 : No Barcode found
        5 : Multiple barcode found
    This list are when we have two good contruct
        6 : Barcorde not the same
        7 : insertitude in multiple barcode on one read after validation
        8 : min one read have a recode error 2 and the rest can't validate a group of barcode
        9 : The rest
    For the test of barcode:
        10 : Multiple match BCX
        11 : Multiple match BCY
        12 : Multiple match BCX and BCY
    Extrait aussi les barcode ainsi que les bonnes sequences"""