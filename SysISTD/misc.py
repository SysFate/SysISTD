'''
Created on 17 Nov 2022

@author: franco
'''

import pandas as pd

def extract_barcode_present(matrix_position):
    dictConvertPositionMatrix = {}
    list_barcode_present = {'x':set(),'y':set()}
    df = pd.read_csv(matrix_position,sep="\t",index_col=0)
    for index, row in df.iterrows():
        for nameCol in list(df):
            if isinstance(row[nameCol],str):
                results = row[nameCol][1:].split('y')
                list_barcode_present['x'].add(results[0])
                list_barcode_present['y'].add(results[1])
                dictConvertPositionMatrix[results[0]+'x'+results[1]] = '{}x{}'.format(index,nameCol)
    return dictConvertPositionMatrix,list_barcode_present,len(df.columns)