'''
Created on 19 Nov 2022

@author: Francois STÜDER

Only feature count was tested in this script
'''

import pandas as pd
import sys

from os import path
from subprocess import Popen, PIPE

from SysISTD import tools

cmdQuantification = "featureCounts -T __core__ -F GTF -t __feature__ -s 0 -a __GTF__ -o __out__"
tsvFolder = "./"

def quantification(samDir, listBarcode, feature, gtfFile, dictConvertPositionMatrix):
    """
    Counting : flexible cmd with config file
    Keywords :
    __GTF__ : GTF file path
    __feature__ : feature use
    __out__ : output path + name of count file
    """
    sys.stdout.write(f'### - Quantification - {feature}\n')
    cmd_count = []
    global tsvFolder
    pathTSVFile = path.join(tsvFolder,f'matrix_count_{feature}.tsv')
    global cmdQuantification
    for partCMD in cmdQuantification.split(' '):
        if partCMD == '__GTF__':
            cmd_count.append(gtfFile)
        elif partCMD == '__core__':
            cmd_count.append(str(tools.core))
        elif partCMD == '__feature__':
            cmd_count.append(feature)
        elif partCMD == '__out__':
            cmd_count.append(pathTSVFile)
        else:
            cmd_count.append(partCMD)
    
    dictConvertsamFileBarcode = {}
    for barcode in listBarcode:
        samFile = path.join(samDir, f'{barcode}.sam')
        cmd_count.append(samFile)
        dictConvertsamFileBarcode[samFile] = barcode
    
    quantiRun = Popen(cmd_count, stdout=PIPE, stderr=PIPE)
    stdout, stderr = quantiRun.communicate()
    with open(path.join(tools.tmp,f"quantification_{feature}.log"), 'w') as log, open(path.join(tools.tmp,f"quantification_{feature}.err"), 'w') as error:
        log.write(stdout.decode('utf-8'))
        error.write(stderr.decode('utf-8'))
    if quantiRun.returncode != 0:
        raise Exception(f'Error with counting {feature}')
    
    countMatrix = pd.read_csv(pathTSVFile, sep = '\t', skiprows=1, index_col=0) # <= skiprows is for the first row of feature count
    summaryCount = pd.read_csv(pathTSVFile+".summary", sep = '\t', index_col=0, header=0)
    #### Clean special feature count ####
    for col in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']:
        try:
            countMatrix.drop(col, axis=1, inplace=True)
        except:
            pass
    countMatrix.rename(dictConvertsamFileBarcode, axis=1, inplace=True)
    countMatrix.rename(dictConvertPositionMatrix, axis=1, inplace=True)
    totalCount = int(countMatrix.to_numpy().sum())
    countMatrix.to_csv(pathTSVFile, sep='\t')
    summaryCount.rename(dictConvertsamFileBarcode, axis=1, inplace=True)
    summaryCount.rename(dictConvertPositionMatrix, axis=1, inplace=True)
    summaryCount.to_csv(pathTSVFile+".summary", sep='\t')
    
    return totalCount