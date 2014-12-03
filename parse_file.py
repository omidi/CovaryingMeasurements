#! /software/bin/python

import csv 
from sys import argv 
import tabix as tb
import numpy as np
import pysam
from dependency_model import dependencyScore


def convert2Matrix(row):
    if len(row)==33:
        tmp = map(float, row[4:])
        return np.array([tmp[0:7], tmp[8:15], tmp[15:22], tmp[22:]])
    elif len(row)==32:
        tmp = map(float, row[3:])
        return np.array([tmp[0:7], tmp[8:15], tmp[15:22], tmp[22:]])
    else:
        return -1
        

def loadRefRegions(infile):
    refRegions = {}
    with open(infile) as inf:
        for record in csv.reader(inf, delimiter="\t"):
            id = '_'.join(record[0:3])
            refRegions.setdefault( id, [] )
            refRegions[ id ] = convert2Matrix(record)
    return refRegions
    

def main():    
    sigma = np.array([0.2, 8., 2., 4.2])
    C = sigma*np.identity(4)
    S = sigma*np.identity(4)    
    offset = 10000  # search within 10 Kb from up/downstream
    if not len(argv) == 3:
        print "Usage: ./prog ref_regions distal_regions.tabix"
        exit()
    refRegions = loadRefRegions(argv[1])
    queryFile = pysam.Tabixfile(argv[2])
    for refRegionId in refRegions.keys():
        X = refRegions[refRegionId]
        chrom, start, end = refRegionId.split('_')
        for distalRegion in queryFile.fetch('chr7', int(start)-offset, int(end)+offset):
            Y = convert2Matrix(distalRegion.split())
            print distalRegion
            print dependencyScore(X, Y, C, S)
            
            
if __name__ == '__main__':
    main()
