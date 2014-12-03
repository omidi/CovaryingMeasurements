#! /software/bin/python

import csv
import sys
import os
import numpy as np


def dotProduct(data, mark1, mark2):
    for v in data['ZT02'].values():
        normalization = 1./(len(v)*7. - 1)
        break
    # pseudo_count = 1.e-6
    result = 0.
    for zeitGeber in ['ZT02', 'ZT06', 'ZT10', 'ZT14', 'ZT18', 'ZT22', 'ZT26']:    
        x_tmp = data[zeitGeber][mark1][ (data[zeitGeber][mark1]>0) & (data[zeitGeber][mark2]>0)]
        x = np.log( x_tmp / np.sum(x_tmp) )
        x_bar = np.mean(x)
        y_tmp = data[zeitGeber][mark2][ (data[zeitGeber][mark1]>0) & (data[zeitGeber][mark2]>0)]
        y = np.log(y_tmp / np.sum(y_tmp))
        y_bar = np.mean(y)
        result += np.dot((x - x_bar), (y - y_bar))
    return normalization*result 
    

def main():
    if len(sys.argv)!=4:
        print "Usage: prog.py data_file mark1 mark2"
        exit(1)
    mark1 = sys.argv[2]
    mark2 = sys.argv[3]
    data = {}
    with open(sys.argv[1]) as inf:        
        for record in csv.DictReader(inf, delimiter='\t'):
            for ZT in ['ZT02', 'ZT06', 'ZT10', 'ZT14', 'ZT18', 'ZT22']:
                rep1_mark = mark1 + ZT
                rep2_mark = mark2 + ZT
                data.setdefault(ZT, {mark1:[], mark2:[]})                
                data[ZT][mark1].append(float(record[rep1_mark]))
                data[ZT][mark2].append(float(record[rep2_mark]))
    dotProd = dotProduct(data, mark1, mark2)
    print dotProd
    

if __name__ == '__main__':
    main()
