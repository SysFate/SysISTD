'''
Created on 17 Nov 2022

@author: Francois STÜDER
'''

import sys

from io import TextIOWrapper,StringIO
from os import path,remove
from subprocess import Popen, PIPE
from threading import Thread,RLock
from time import time

from SysISTD import tools
from SysISTD import alignment
from SysISTD import AnnalyseRecord
from SysISTD import quantification

ID = None
cmdDecompress = ["pigz", "-dc"]
errorGuib = 1
errorBC = 1
calculPool = None
genomePath = ""
gtfFile = ""
featuresToGenerate = ["exon","transcript"]
posMatrix = ""
sizeX = 0

def demultiple(fastq_file_1, fastq_file_2, outFolder, dictConvertPositionMatrix, verbose):
    """
    Demultiplexing : get coordinates (X & Y) from fastq file
    - fastq_file_1     : fastq file 1
    - fastq_file_2     : fastq file 2
    - sequence             : reference sequence for get X and Y sequences (In this case Gibson sequence)
    - sequence_reverse    : reverse reference sequence for get reverse X and reverse Y sequences (In this case reverse Gibson sequence)
    """
    sys.stdout.write('### - Demultiplexing\n')
    statsOut = {-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,'total_read_count':0,'goodPolyT':0,'goodG':0,'reads_good_bacode':0}
    dicCountFindReadByCoordinate = {}
    UMIStats = {'number_by_barcode':{},'number_by_pos':{},'mean_by_barcode':{}}
    
    if fastq_file_1[-2:] == 'gz' and fastq_file_2[-2:] == 'gz':
        #with gzip.open(fastq_file_1, 'rt') as read1, gzip.open(fastq_file_2, 'rt') as read2:
        #    searchTarget(read1, read2, dic_find, list_test, dic_X, dic_Y)
        global cmdDecompress
        cmdreadFile1 = cmdDecompress.copy()
        cmdreadFile1.append(fastq_file_1)
        cmdreadFile2 = cmdDecompress.copy()
        cmdreadFile2.append(fastq_file_2)
        readFile1 = Popen(cmdreadFile1, stdout=PIPE, stderr=PIPE)
        readFile2 = Popen(cmdreadFile2, stdout=PIPE, stderr=PIPE)
        cmdreadFile1 = None
        cmdreadFile2 = None
        clean_arrange_fastq_by_coordinate(TextIOWrapper(readFile1.stdout, encoding="utf-8"), TextIOWrapper(
            readFile2.stdout, encoding="utf-8"), statsOut, dicCountFindReadByCoordinate, UMIStats, verbose)
    else:
        with open(fastq_file_1, 'rt') as read1, open(fastq_file_2, 'rt') as read2:
            clean_arrange_fastq_by_coordinate(read1, read2, statsOut, dicCountFindReadByCoordinate, UMIStats, verbose)
            
    total_read_count = statsOut['total_read_count'] #total number of reads
    reads_good_bacode = statsOut['reads_good_bacode'] #reads with structure of probe and good barcode(s)
    goodPolyT = statsOut['goodPolyT'] + statsOut['reads_good_bacode']
    goodG = statsOut['goodG'] + goodPolyT
    
    sys.stdout.write("\n"+"#"*25+"\n")
    sys.stdout.write("#"*7+" A rapid overview of all stats "+"#"*7+"\n")
    sys.stdout.write("Total number of reads {}\n".format(total_read_count))
    sys.stdout.write("Total Gibson : {} | {} %\n".format(goodG, (goodG/total_read_count)*100))
    sys.stdout.write("Good construct {}\n".format(goodPolyT))
    sys.stdout.write("Number of coordinate found {}\n".format(len(dicCountFindReadByCoordinate)))
    sys.stdout.write("#"*25+"\n")
    if verbose:
        sys.stdout.write("\n"+"#"*25+"\n")
        sys.stdout.write("#"*7+" Number reads by coordinate "+"#"*7+"\n")
        sys.stdout.write(str(dicCountFindReadByCoordinate))
        sys.stdout.write("\n"+"#"*25+"\n")
        
        sys.stdout.write("\n"+"#"*25+"\n")
        sys.stdout.write("#"*7+" Code errors of sequences "+"#"*7)
        sys.stdout.write(AnnalyseRecord.get_stats_code_sequences())
        sys.stdout.write("\n"+"="*25+"\n")
        sys.stdout.write(str(statsOut))
        sys.stdout.write("\n"+"#"*25+"\n")
        sys.stdout.write("#"*7+" UMIStats V2 "+"#"*7+"\n")
        for nameUMIStat, statsUMIToPrint in UMIStats.items():
            sys.stdout.write(str(nameUMIStat)+":\n")
            sys.stdout.write(str(statsUMIToPrint))
            sys.stdout.write("\n"+"="*25+"\n")
        sys.stdout.write("#"*25+"\n")
    UMIStats = None
    statsOut = None
    
    readMapped = generate_sam()
    featureCountStats = generation_quantification_heatmap(dicCountFindReadByCoordinate.keys(),dictConvertPositionMatrix,outFolder)
    
    logfileWriter(total_read_count, goodG, goodPolyT, reads_good_bacode, readMapped, featureCountStats, dicCountFindReadByCoordinate, outFolder)
    

def clean_arrange_fastq_by_coordinate(read1, read2, statsOut, dicCountFindReadByCoordinate, UMIStats, verbose):
    sys.stdout.write('### - Clean and arrange fastq\n')
    
    UMIStatsByBarcode = {}
    # On ouvre les fichiers pour les parcourir
    read1Parse = read1#read1Parse = SeqIO.parse(read1,"fastq")
    read2Parse = read2#read2Parse = SeqIO.parse(read2,"fastq")
    # On va garder en memoire les reads qui n'ont pas le bon copain
    read1Alone = {}
    read2Alone = {}
    # Cette etape est parallelise pour plus de vitesse
    #listReadsByJob = []
    listThreads = tools.init_list_thread(tools.core)
    listLaunchThreads = tools.init_list_thread(int(int(tools.core)/1.5)+1)
    # On se stop lorsque record1 ou record2 est au bout
    recordReadable = True
    record1Readable = True
    record2Readable = True
    record1Lecture = extract_record_file(read1Parse)
    record1Lecture.start()
    record2Lecture = extract_record_file(read2Parse)
    record2Lecture.start()
    statTimeTake = {'get_reads':0,'all_the_loop_before_add_in_queue':0,'all_the_loop':0,
                    'wait_reads':0,'before_lock_reads_alone':0,'before_lock_launch_annalyse':0,'after_launch_annalyse':0,
                    'get_listRead_lock':0,'get_stats_lock':0,'get_UMIStats_lock':0,'get_queue_lock':0,'all_analyse':0,
                    'numberItteration':0}
    while recordReadable:
        beginWait = time()
        statTimeTake['numberItteration'] += 1
        record1Lecture.join()
        record2Lecture.join()
        statTimeTake['get_reads'] += (time()-beginWait)
        dicoRecord1 = record1Lecture.groupLines #record1Lecture.outDict
        dicoRecord2 = record2Lecture.groupLines #record2Lecture.outDict
        if record1Lecture.alreadyReadable and record2Lecture.alreadyReadable:
            record1Lecture = extract_record_file(read1Parse)
            record1Lecture.start()
            record2Lecture = extract_record_file(read2Parse)
            record2Lecture.start()
        else:
            record1Readable = record1Lecture.alreadyReadable
            record2Readable = record2Lecture.alreadyReadable
            recordReadable = False
        
        statTimeTake['all_the_loop_before_add_in_queue'] += (time()-beginWait)
        tools.auto_add_thread(listLaunchThreads, 0.0001,Thread(target=prepare_and_launch_reads_treatment, args=(dicoRecord1, dicoRecord2, read1Alone, read2Alone, listThreads, dicCountFindReadByCoordinate, statsOut, UMIStatsByBarcode,statTimeTake)))
        statTimeTake['all_the_loop'] += (time()-beginWait)

    dicoRecord1 = None
    dicoRecord2 = None
    tools.wait_end_thread_list(listLaunchThreads)

    numberItteration = statTimeTake['numberItteration']
    del statTimeTake['numberItteration']
    statTimeTake = {key: value / numberItteration for key, value in statTimeTake.items()}
    sys.stdout.write(str(statTimeTake))

    listReadsByJob = []

    try:
        while record1Readable:
            record = next(read1Parse)
            recordID = record.id.split('/')[0]
            if recordID in read2Alone:
                listReadsByJob.append([record,read2Alone[recordID]])
                del read2Alone[recordID]
            else:
                read1Alone[recordID] = record
    except:
        pass
    try:
        while record2Readable:
            record = next(read2Parse)
            recordID = record.id.split('/')[0]
            if recordID in read1Alone:
                listReadsByJob.append([read1Alone[recordID],record])
                del read1Alone[recordID]
            else:
                read2Alone[recordID] = record
    except:
        pass
    
    if len(listReadsByJob) != 0:
        tools.auto_add_thread(listThreads,0.001,AnnalyseRecord.AnnalyseRecord(listReadsByJob,dicCountFindReadByCoordinate,statsOut,UMIStatsByBarcode))
    listReadsByJob = None
    record = None
        
    sys.stdout.write("\n#########################\nRead alone 1: "+str(len(read1Alone))+"\nRead alone 2: "+str(len(read2Alone))+"\n#########################\n")
    read1Alone = None
    read2Alone = None
    
    tools.wait_end_thread_list(listThreads)

    mergedDicFastqCoordinate = {}
    AnnalyseRecord.merge_dic_fastq_coordinate(mergedDicFastqCoordinate, AnnalyseRecord.queueWriteGoodRead)
    AnnalyseRecord.queueWriteGoodRead = None
    AnnalyseRecord.write_coordinate_fastq(mergedDicFastqCoordinate)
    mergedDicFastqCoordinate = None
    
    ##################################################################################################
    ######################################## WORK In progress ########################################
    for barcode,data in UMIStatsByBarcode.items():
        UMIStats['mean_by_barcode'][barcode] = {}
        UMIStats['number_by_barcode'][barcode] = {}
        for pos,data2 in data.items():
            # On a besoin de la moyenne par barcode
            numberUMIBCHere = []
            for seq,value in data2.items():
                numberUMIBCHere.append(value)
            UMIStats['mean_by_barcode'][barcode][pos] = sum(numberUMIBCHere)/len(numberUMIBCHere)
            if pos in UMIStats['number_by_pos']:
                UMIStats['number_by_pos'][pos].append(len(data2))
            else:
                UMIStats['number_by_pos'][pos] = [len(data2)]
            UMIStats['number_by_barcode'][barcode][pos] = len(data2)
    if verbose:
        sys.stdout.write("\n"+"#"*25+"\n"+"#"*7+" UMI Stats by barcode "+"#"*7+"\n")
        sys.stdout.write(str(UMIStatsByBarcode))
        sys.stdout.write("\n"+"#"*25+"\n")
    ##################################################################################################

class extract_record_file(Thread):
    '''
    classdocs
    '''


    def __init__(self, fileHandleRecord):
        '''
        Constructor
        '''
        Thread.__init__(self)
        self.fileHandleRecord = fileHandleRecord
        self.groupLines = None
        self.alreadyReadable = True

    def run(self):
        self.groupLines = []
        try:
            for h in range(0, int((tools.core-1)/4)+1):
                lines = []
                self.groupLines.append(lines)
                for i in range(0, AnnalyseRecord.readPas*16):
                    lines.append(next(self.fileHandleRecord))
        except:
            if len(self.groupLines[-1]) == 0:
                self.groupLines = self.groupLines[:-1]
            self.alreadyReadable = False
            
verrouReadAlone = RLock()
verrouLaunchAnnalyse = RLock()
def prepare_and_launch_reads_treatment(dicoRecord1, dicoRecord2, read1Alone, read2Alone, listThreads, dicCountFindReadByCoordinate, statsOut, UMIStatsByBarcode, statTimeTake):
    beginWait = time()
    record1Convert = convert_list_record_file(dicoRecord1)
    record1Convert.start()
    record2Convert = convert_list_record_file(dicoRecord2)
    record2Convert.start()
    record1Convert.join()
    dicoRecord1 = record1Convert.dicoRecord
    record2Convert.join()
    dicoRecord2 = record2Convert.dicoRecord
    statTimeTake['wait_reads'] += (time()-beginWait)
    
    listReadsByJob = []
    potentialReadAlone1 = {}
    for recordID, record in dicoRecord1.items():
        if recordID in dicoRecord2:
            listReadsByJob.append([record, dicoRecord2[recordID]])
            del dicoRecord2[recordID]
        else:
            potentialReadAlone1[recordID] = record
    
    record1Convert = None
    record2Convert = None
    statTimeTake['before_lock_reads_alone'] += (time()-beginWait)
    with verrouReadAlone:
        for recordID, record in potentialReadAlone1:
            if recordID in read2Alone:
                listReadsByJob.append([record, read2Alone[recordID]])
                del read2Alone[recordID]
            else:
                read1Alone[recordID] = record
    
        for recordID, record in dicoRecord2.items():
            if recordID in read1Alone:
                listReadsByJob.append([read1Alone[recordID], record])
                del read1Alone[recordID]
            else:
                read2Alone[recordID] = record
    dicoRecord1 = None
    dicoRecord2 = None
    potentialReadAlone1 = None
    statTimeTake['before_lock_launch_annalyse'] += (time()-beginWait)
    with verrouLaunchAnnalyse:
        tools.auto_add_thread(listThreads, 0.0001, AnnalyseRecord.AnnalyseRecord(listReadsByJob, dicCountFindReadByCoordinate, statsOut, UMIStatsByBarcode, statTimeTake))
    statTimeTake['after_launch_annalyse'] += (time()-beginWait)
    
class convert_list_record_file(Thread):
    '''
    classdocs
    '''


    def __init__(self, groupLines):
        '''
        Constructor
        '''
        Thread.__init__(self)
        self.groupLines = groupLines
        self.dicoRecord = None

    def run(self):
        global calculPool
        jobsRecordConvert = []
        for lines in self.groupLines:
            jobsRecordConvert.append(calculPool.apply_async(convert_record_file, (StringIO("".join(lines)),)))
        self.dicoRecord = {}
        for recordConvert in jobsRecordConvert:
            self.dicoRecord.update(recordConvert.get())

def convert_record_file(fileHandleRecord):
    from Bio import SeqIO
    outDict = {}
    try:
        seqIOHandler = SeqIO.parse(fileHandleRecord,"fastq")
        while True:
            record = next(seqIOHandler)
            outDict[record.id.split('/')[0]] = record
    except:
        pass
    return outDict

def generate_sam():
    # Merge all little fastq maybe useless since I add the rename on each sequence
    goodSeqFastqFile1 = path.join(tools.tmp,"validate_seq_1.fastq")
    mergeFastqFile1 = Popen("cat {} > {}".format(path.join(AnnalyseRecord.fastqDir,"*_1.fastq"),goodSeqFastqFile1),shell=True, stdout=PIPE, stderr=PIPE)
    goodSeqFastqFile2 = path.join(tools.tmp,"validate_seq_2.fastq")
    mergeFastqFile2 = Popen("cat {} > {}".format(path.join(AnnalyseRecord.fastqDir,"*_2.fastq"),goodSeqFastqFile2),shell=True, stdout=PIPE, stderr=PIPE)
    # Prepare some stats ?
    global genomePath
    
    # Validate the merge
    stdout, stderror = mergeFastqFile1.communicate()
    if mergeFastqFile1.returncode != 0:
        raise Exception("Error during the fastq merge: {}".format(stderror))
    stdout, stderror =mergeFastqFile2.communicate()
    if mergeFastqFile2.returncode != 0:
        raise Exception("Error during the fastq merge: {}".format(stderror))
    # Launch the alignment and fastq generation
    sys.stdout.write('\n### - Alignment + Generate sam\n')
    readMapped = alignment.alignment_generate_sam_coordinate(goodSeqFastqFile1, goodSeqFastqFile2, genomePath)
    # Remove the two merged fastq
    remove(goodSeqFastqFile1)
    remove(goodSeqFastqFile2)
    
    return readMapped
    
def generation_quantification_heatmap(listBarcode,dictConvertPositionMatrix,outFolder):
    global gtfFile
    global featuresToGenerate
    featureCountStats = {}
    for feature in featuresToGenerate:
        featureCountStats[feature] = quantification.quantification(alignment.samDir, listBarcode, feature, gtfFile, dictConvertPositionMatrix)
        try:
            heatmapFromMatrix(path.join(quantification.tsvFolder,f'matrix_count_{feature}.tsv'), path.join(outFolder,f'matrix_heatmap_{feature}.tsv'))
        except Exception as e:
            pass
    return featureCountStats

def heatmapFromMatrix(matrixCount, matrixHeatMap):
    """
    matrixCount : matrix count filter
    matrix_position : all position (coordiante) in the correct order on the glass slide
    matrixHeatMap : output matrix for heatmap visualization
    """
    global posMatrix
    global sizeX
    with open(posMatrix, 'r') as position, open(matrixCount, 'r') as count, open(matrixHeatMap, 'w') as output:
        dicIndex = {}
        dic = {}
        index_line = 0
        for line in count:
            if index_line == 0:
                for i in range(len(line[:-1].split('\t'))):
                    key = line[:-1].split('\t')[i]
                    if len(key) > 0:
                        keyFormat1 = key.split('x')[0]
                        keyFormat2 = key.split('x')[1]
                        dicIndex[i] = f'x{keyFormat1}y{keyFormat2}'
            else:
                for i in range(len(line[:-1].split('\t'))):
                    if i != 0:
                        if int(line[:-1].split('\t')[i]) > 0:
                            if dicIndex[i] not in dic.keys():
                                dic[dicIndex[i]] = 1
                            else:
                                dic[dicIndex[i]] = dic[dicIndex[i]]+1
            index_line += 1
        index_row = 0
        indexFitMap = 0
        for line in position:
            if not line.split('\t')[0].isdigit():
                tmp_index = ''
                for o in range(sizeX):
                    tmp_index = '{}\t{}'.format(tmp_index, o+1)
                output.write('{}'.format(tmp_index))
            else:
                output.write('{}'.format(index_row))
                for o in range(sizeX):
                    if line.split('\t')[o+1] in dic.keys():
                        output.write('\t{}'.format(dic[line[:-1].split('\t')[o+1]]))
                        indexFitMap += 1
                    else:
                        output.write('\tNA')
            output.write('\n')
            index_row += 1
        output.close()

def logfileWriter(total_read_count, reads_target_sequence, reads_structure_probe, reads_good_bacode, reads_mapped, featureCountStats, find_coordinates, outFolder):
    """
    logfileWriter gives informations about reads in logfile (text format)
    Which number of reads on each step of the pipeline
    total_read_count :                 The total number of read in fastq file
    reads_target_sequence :            Reads which have the target sequence
    reads_structure_probe :            Reads which have the correct structure of the probe (extact distance)
    reads_good_bacode :                Reads which have the good structure and barcodes match with the input barcodes' list
    reads_mapped :                    Reads which are mapped on a reference genome
    """
    global ID
    with open(path.join(outFolder,"demultiplexing.log"), 'w') as logfile:
        logfile.write(f'####################################\n#\t\t\tLogfile - {ID}\n####################################\n\n')
        logfile.write(f'number of read : {total_read_count}\n')
        logfile.write(f'STEP 1 - number of reads with target sequence : {reads_target_sequence}\t| {(reads_target_sequence*100)/total_read_count} %\n')
        logfile.write(f'STEP 2 - number of reads with probe structure : {reads_structure_probe}\t| {(reads_structure_probe*100)/total_read_count} %\n')
        logfile.write(f'STEP 3 - number of reads with probe structure and good barcodes : {reads_good_bacode}\t| {(reads_good_bacode*100)/total_read_count} %\n')
        logfile.write(f'STEP 4 - number of reads mapped : {int(reads_mapped)}\t| {(int(reads_mapped)*100)/total_read_count} %\n')
        numberPrime = 0
        mapFeature = []
        for feature, count in featureCountStats.items():
            prime = "'"*numberPrime
            logfile.write(f'STEP 5{prime} - number of reads count {feature} : {count}\t| {((count*100)/total_read_count)} %\n')
            logfile.write(f'\tIt remains {int(reads_mapped-count)} which are mapped on genome but are not in {feature}\n')
            mapFeature.append([feature, count])
            numberPrime += 1
        for i in range(len(mapFeature)):
            for j in range(i+1,len(mapFeature)):
                logfile.write(f'\nWe have {mapFeature[i][1]-mapFeature[j][1]} reads are in {mapFeature[j][0]} but not in {mapFeature[i][0]}\n')
        logfile.write(f'\nNumber of gexel find in the demultiplexing step : {len(find_coordinates)}\n')