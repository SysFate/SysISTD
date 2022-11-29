'''
Created on 24 Apr 2020

@author: Francois STÜDER

'''

import argparse
import sys

from Bio.Seq import Seq
from datetime import datetime
from multiprocessing import Pool
from os import path
from time import time

from SysISTD import __version__
from SysISTD import tools
from SysISTD import misc
from SysISTD import alignment
from SysISTD import AnnalyseRecord
from SysISTD import quantification
from SysISTD import demultiplex

def create_coordinates():
    """
    Create dictionnary from a tsv file (arg3).
    Get X and Y from tsv file and generate reverse sequences of X and Y.
    - coordinates     : tsv file all X and Y sequences
    TODO : use coordinateObject
    """
    line_count = 1
    inFile = open(coordinates, "r")
    for line in inFile:
        cutLine = line.split('\t')
        if cutLine[0].isdigit():
            value = cutLine[0]
            keyX = cutLine[1]
            keyY = cutLine[2][:-1]
            line_count += 1
            if len(keyX) > AnnalyseRecord.lengthBC:
                AnnalyseRecord.lengthBC = len(keyX)
            if len(keyY) > AnnalyseRecord.lengthBC:
                AnnalyseRecord.lengthBC = len(keyY)
            if value in list_barcode_present['x']:
                AnnalyseRecord.barcodesX[value] = {'+' : '('+keyX+'){s<='+str(demultiplex.errorBC)+'}', '-' : '('+str(Seq(keyX).reverse_complement())+'){s<='+str(demultiplex.errorBC)+'}'}
            if value in list_barcode_present['y']:
                AnnalyseRecord.barcodesY[value] = {'+' : '('+str(Seq(keyY).reverse_complement())+'){s<='+str(demultiplex.errorBC)+'}', '-' : '('+keyY+'){s<='+str(demultiplex.errorBC)+'}'}

def main(argv=None):
    begin = time()
    parser = argparse.ArgumentParser(description=f'This software developed by the SysFate team process FASTQ files by extracting sequences from each barcode position.\n You have the {__version__}', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument("--f1", metavar='file1', type=str, required=True, help="The file 1 you want analyse")
    parser.add_argument("--f2", metavar='file2', type=str, required=True, help="The file 2 you want analyse")
    parser.add_argument("-n","--name", metavar='name', type=str, default=datetime.now().strftime("%Y-%m-%d_%H:%M:%S"), help="You can give a name for your analyse. By default the software use the begin date.")
    parser.add_argument("--tmp", metavar='tmp', type=str, default="/tmp", help="Folder where send temporary files")
    parser.add_argument("-o","--out", metavar='outdir', type=str, default=".", help="Folder where send files")
    parser.add_argument("--coor", metavar='coordinate', type=str, required=True, help="TSV file who contain coordinates code with nucleique barcode")
    parser.add_argument("--pos", metavar='position', type=str, required=True, help="TSV file who contain coordinates code with their position in the grid")
    parser.add_argument("--guibson", metavar='guibson', type=str, default="ACATTGAAGAACCTGTAGATAACTCGCTGT", help="Sequences guibson")
    parser.add_argument("-G","--guibsonE", metavar='guibsonE', type=int, default=1, help="Guibson error accepted")
    parser.add_argument("-B","--barcodeE", metavar='barcodeE', type=int, default=1, help="Barcode error accepted")
    parser.add_argument("--polyTE", metavar='polyTE', type=int, default=1, help="PolyT error accepted")
    parser.add_argument("-c","--core", metavar='core', type=int, default=1, help="Number of core the software can use")
    parser.add_argument("--polyTSizeMin", metavar='polyTSizeMin', type=int, default=12, help="Size min for a polyT")
    parser.add_argument("--nopolytsearch", dest='nopolytsearch', default=False, action='store_true', help="If you don't want validate the research with polyT")
    parser.add_argument("--trimSequence", metavar='trimSequence', type=str, default="", help="After the UMI, what do you want trim ? A fix number of base or a specific sequence. By default, trim the polyT as you config before")
    parser.add_argument("--noUMIsearch", dest='noUMIsearch', default=False, action='store_true', help="If you don't want research UMI")
    parser.add_argument("-a","--aligner", metavar='aligner', type=str, default="bowtie2 --reorder --very-sensitive -p __core__ -x __genome__ -1 __file1__ -2 __file2__ -S __outputSam__", help="This is the aligner command. You can write what you want if the __file1__, __file2__, __genome__ is provide and the output is configured to send the results in __outputSam__.")
    parser.add_argument("-g","--genome", metavar='genome', type=str, required=True, help="Aligner genome index path")
    parser.add_argument("-q","--quantifier", metavar='quantifier', type=str, default="featureCounts -T __core__ -F GTF -t __feature__ -s 0 -a __GTF__ -o __out__", help="This is the read quantifier command. You can write what you want if the __GTF__, __feature__ is provide and the output is configured to send the results in __out__.")
    parser.add_argument("--feature", metavar='feature', type=str, default="exon,transcript", help="Feature use during the read quantifier in your file separate by a comma")
    parser.add_argument("--gtf", metavar='gtf', type=str, required=True, help="GTF file to use with feature count")
    parser.add_argument("--decompress", metavar='decompress', type=str, default="pigz -dc", help="Command use to decompress fastq files")
    parser.add_argument("--dev", dest='dev', default=False, action='store_true', help="Keep all temporar datas")
    args = parser.parse_args()
    if args.f1 == args.f2:
        raise Exception(f"You gave in file 1 and file 2 the same file !")
    tools.test_file_can_open(args.f1, "Fastq file 1 file can't be open")
    tools.test_file_can_open(args.f2, "Fastq file 2 file can't be open")
    tools.test_file_can_open(args.gtf, "GTF file can't be open")
    global coordinates
    coordinates = args.coor

    demultiplex.ID = args.name
    outFolder = path.join(args.out,args.name)
    if path.exists(outFolder):
        i = 0
        while path.exists(outFolder+"."+str(i)):
            i = i+1
        outFolder = outFolder+"."+str(i)
    tools.make_dirs(outFolder)
    sys.stdout.write(f'### - {args.name}\n')
    sys.stdout.write(f'Sotware version - {__version__}\n')
    tools.tmp = path.join(args.tmp, f'{args.name}_demultiplex_tmp')
    if path.exists(tools.tmp):
        i = 0
        while path.exists(tools.tmp+"."+str(i)):
            i = i+1
        tools.tmp = tools.tmp+"."+str(i)
    tools.make_dirs(tools.tmp)
    AnnalyseRecord.fastqDir = path.join(tools.tmp, 'fastq_files')
    AnnalyseRecord.outBadSeqFolder = path.join(tools.tmp,f"fastq_bad_sequences")
    tools.make_dirs(AnnalyseRecord.fastqDir)
    tools.make_dirs(AnnalyseRecord.outBadSeqFolder)
    
    demultiplex.cmdDecompress = args.decompress.split(" ")
    demultiplex.errorGuib = args.guibsonE
    demultiplex.errorBC = args.barcodeE
    tools.core = args.core
    AnnalyseRecord.guibsonSeq = {'+' : '(?e)('+args.guibson+'){s<='+str(demultiplex.errorGuib)+'}', '-' : '(?e)('+str(Seq(args.guibson).reverse_complement())+'){s<='+str(demultiplex.errorGuib)+'}'}

    # Configuration for polyT
    AnnalyseRecord.polyTSearch = not args.nopolytsearch
    ##################################################################################################
    ############ This Part will not be use by the main analyser if polyTSearch is disabled ########### 
    if args.polyTE > 0:
        polySearch = {"A": '(?e)(A{'+str(args.polyTSizeMin)+'}A+){s<='+str(args.polyTE)+'}', "T": '(?e)(T{'+str(args.polyTSizeMin)+'}T+){s<='+str(args.polyTE)+'}'}
    else:
        polySearch = {"A": 'A{'+str(args.polyTSizeMin)+'}A+', "T": 'T{'+str(args.polyTSizeMin)+'}T+'}

    if args.polyTSizeMin > 0:
        AnnalyseRecord.polyT = {"+": AnnalyseRecord.guibsonSeq['+']+'.{14,20}'+polySearch['T'], "-": polySearch['A']+'.{14,20}'+AnnalyseRecord.guibsonSeq['-'], "s+": polySearch['T']+'.{14,20}'+AnnalyseRecord.guibsonSeq['+'], "s-": AnnalyseRecord.guibsonSeq['-']+'.{14,20}'+polySearch['A']}
    else:
        AnnalyseRecord.polyT = {"+": AnnalyseRecord.guibsonSeq['+']+'.{14,20}', "-": '.{14,20}'+AnnalyseRecord.guibsonSeq['-'], "s+": '.{14,20}'+AnnalyseRecord.guibsonSeq['+'], "s-": AnnalyseRecord.guibsonSeq['-']+'.{14,20}'}
    ##################################################################################################

    ##################################################################################################
    ######################################## WORK In progress ########################################
    if args.trimSequence.isdigit():
        AnnalyseRecord.trimSequence = int(args.trimSequence)
        AnnalyseRecord.polyTPair = {"+": polySearch['T'], "-": polySearch['A']}
    elif isinstance(args.trimSequence,str) and len(args.trimSequence) > 0:
        AnnalyseRecord.trimSequence = args.trimSequence
        AnnalyseRecord.polyTPair = {"+": polySearch['T'], "-": polySearch['A']}
    else:
        AnnalyseRecord.trimSequence = 0
        AnnalyseRecord.polyTPair = {"+": polySearch['T'], "-": polySearch['A']}
        global UMISearch
    AnnalyseRecord.UMISearch = not args.noUMIsearch
    ##################################################################################################
    polySearch = None
    
    demultiplex.genomePath = args.genome
    alignment.cmdAlignment = args.aligner
    alignment.samDir = path.join(tools.tmp, 'sam_files')
    tools.make_dirs(alignment.samDir)
    
    demultiplex.gtfFile = args.gtf
    demultiplex.featuresToGenerate = args.feature.split(",")
    quantification.cmdQuantification = args.quantifier
    quantification.tsvFolder = outFolder
    
    demultiplex.posMatrix = args.pos
    
    global list_barcode_present
    dictConvertPositionMatrix,list_barcode_present,sizeX = misc.extract_barcode_present(args.pos)
    create_coordinates()
    demultiplex.sizeX = sizeX
    AnnalyseRecord.calculPool = Pool(processes=tools.core-1)
    demultiplex.calculPool = AnnalyseRecord.calculPool
    list_barcode_present = None
    demultiplex.demultiple(args.f1, args.f2, outFolder, dictConvertPositionMatrix, args.dev)
    
    if args.dev:
        tools.move_dir(tools.tmp, outFolder)
    else:
        tools.remove_dir(tools.tmp)
    
    sys.stdout.write("\n The software take {} to be execute".format(time()-begin))

if __name__ == '__main__':
    main()