#! /software/bin/python

import csv
import sys
import os
import argparse
import pysam
import numpy as np
from fitting_codependencies import dotProduct
from itertools import combinations_with_replacement
from fitting_variance import mleFit, empiricalFit
from dependency_model import dependencyScore, flatLikelihoodScore


TIMES=['ZT02', 'ZT06', 'ZT10', 'ZT14', 'ZT18', 'ZT22', 'ZT26']
MARKS=['Pol2', 'H3K4me3', 'DNase', 'H3K27ac']
INDEX={'Pol2':0, 'H3K4me3':1, 'DNase':2, 'H3K27ac':3}

def arguments():
    parser = argparse.ArgumentParser(description="""
    Identifies statistical dependency between genomic regions,
    of which a subset of them labled as Proximal and TSS and
    the rest as Distal.
    For each region, there are N measurements over T time points.
    The arguments are:
    an input filename that contains at each line the location of the signal
    and also the measurements.
    Another file is a tabix indexed file for ONLY distal or enhancer regions. 
    """)
    parser.add_argument('-i', '--datafile', dest='datafile', action='store',
                       type=str, required=True,
                       help='A file that contains the total data.')
    parser.add_argument('-a', '--anchor', dest='anchorfile', action='store',
                       type=str, required=True,
                       help='A file that contains the anchor regions.')    
    parser.add_argument('-d', '--distal', dest='distalfile', action='store',
                       type=str, required=True,
                       help='A file that contains distal file in .gz format (tabix).')
    parser.add_argument('-g', '--header', dest='headerfile', action='store',
                       type=str, required=True,
                       help='The header of the data file.')    
    args = parser.parse_args()
    return args



def loadData(datafile):
    data = {}
    with open(datafile) as inf:        
        for record in csv.DictReader(inf, delimiter='\t'):
            for ZT in TIMES:
                data.setdefault(ZT, {})
                for mark in MARKS:
                    data[ZT].setdefault(mark, [])
                    data[ZT][mark].append(float(record[mark + ZT]))
    for ZT in TIMES:
        for mark in MARKS:
            data[ZT][mark] = np.log2((np.array(data[ZT][mark]) + 1.0))    # plus 1.0 is the pseudo-count
            data[ZT][mark] -= np.mean(data[ZT][mark])
            data[ZT][mark] /= np.std(data[ZT][mark])
    return data



def convertToSlope(X):
    return (X[:, 1:] - X[:, :-1])


def calculateCorrelationMatrix(data):
    correlationMatrix = np.zeros(len(MARKS)**2).reshape(len(MARKS), len(MARKS))
    for marks in combinations_with_replacement(MARKS, 2):
        correlationMatrix[INDEX[marks[0]]][INDEX[marks[1]]] = \
          dotProduct(data, marks[0], marks[1])
        correlationMatrix[ INDEX[marks[1]] ] [ INDEX[marks[0]] ] = \
          correlationMatrix [ INDEX[marks[0]] ] [ INDEX[marks[1]] ]
    return correlationMatrix


def calculateCovarianceMatrix(data):
    covarianceMatrix = np.zeros(len(MARKS)**2).reshape(len(MARKS), len(MARKS))
    for mark in MARKS:
        index = INDEX[mark]
        # rep1, rep2 = [], []
        # for a, b in zip(data[TIMES[0]][mark], data[TIMES[-1]][mark]):
        #     if a > 5. and b > 5.:
        #         rep1.append(a)
        #         rep2.append(b)
        covarianceMatrix[index][index] = \
          empiricalFit( data[TIMES[0]][mark], data[TIMES[-1]][mark] )
    return covarianceMatrix


def loadHeader(headerfile):
    with open(headerfile) as inf:
        header = inf.readline().split()
    return dict([(header[i], i) for i in xrange(len(header))])
    

def gaussian_weight(x):
    return np.exp( -0.5*x**2 )


def convert2MatrixProximal(record, queryFile, width, headerIndex, Mean, SD):
    mat = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))
    records = [record]
    weights = [1.]
    start = int(record[ headerIndex['DHSpeak'] ])
    if (start-2000) < 0:
        start = 2000
    for distalRegion in queryFile.fetch(record[ headerIndex['chrom'] ], start-2500, start+2500):
        record = distalRegion.split()
        if int(record[1])==start:  # if it's excatly the same record
            continue
        records.append(record)
        weight = gaussian_weight(np.abs(start - int(record[1]))/width)
        weights.append(weight)
    weights = np.array(weights) / np.sum(weights)
    for record, weight in zip(records, weights):
        for mark, m in zip(MARKS, range(len(MARKS))):
            for ZT, t in zip(TIMES, range(len(TIMES))):
                mat[INDEX[mark]][t] += weight*float(record[ headerIndex[mark + ZT] ])
    X = np.log2(mat + 1.0)
    X -= Mean
    X /= SD
    return X


def convert2MatrixDistal(row, queryFile, width, headerIndex, Mean, SD):
    mat = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))
    record = row.split()
    records = [record]
    weights = [1.]    
    start = int(record[ headerIndex['DHSpeak'] ])
    if (start-2000) < 0:
        start = 2000
    for distalRegion in queryFile.fetch(record[ headerIndex['chrom'] ], start-2500, start+2500):       
        record = distalRegion.split()
        if int(record[1])==start:  # if it's excatly the same record
            continue
        records.append(record)
        weight = gaussian_weight(np.abs(start - int(record[1]))/width)
        weights.append(weight)        
    weights = np.array(weights) / np.sum(weights)
    for record, weight in zip(records, weights):
        for mark, m in zip(MARKS, range(len(MARKS))):
            for ZT, t in zip(TIMES, range(len(TIMES))):
                mat[INDEX[mark]][t] += weight*float(record[ headerIndex[mark + ZT] ])
    # M = (np.log(mat) - Normalization)
    X = np.log2(mat + 1.0)
    X -= Mean
    X /= SD    
    return X


def calculateNormalization(data):
    normalizationMatrix = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))
    for mark, m in zip(MARKS, range(len(MARKS))):
        for ZT, t in zip(TIMES, range(len(TIMES))):
            normalizationMatrix[m][t] = np.log(np.sum( data[ZT][mark]))
    return normalizationMatrix


def calculateMeanSdMatrix(datafile):
    data = {}
    with open(datafile) as inf:
        for record in csv.DictReader(inf, delimiter='\t'):
            for ZT in TIMES:
                data.setdefault(ZT, {})
                for mark in MARKS:
                    data[ZT].setdefault(mark, [])
                    data[ZT][mark].append(float(record[mark + ZT]))
    M = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))
    S = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))
    for mark, m in zip(MARKS, range(len(MARKS))):
        for ZT, t in zip(TIMES, range(len(TIMES))):
            x = np.log2((np.array(data[ZT][mark]) + 1.0))
            M[m][t] = np.mean(x)
            S[m][t] = np.std(x)
    return M, S


def calculateSdMatrix(data):
    M = np.zeros(len(TIMES)*len(MARKS)).reshape(len(MARKS), len(TIMES))    
    for mark, m in zip(MARKS, range(len(MARKS))):
        for ZT, t in zip(TIMES, range(len(TIMES))):
            v = data[ZT][mark]
            M[m][t] = np.std(v) 
            # M[m][t] = np.sum(v)/len(v)
    return M


def plot(X,Y, name):
    f = open(name + '.R', 'w')
    f.write('\n'.join([
        "d = read.table('%s')" % name,
        "pdf('%s.pdf', height=5, width=7)" % name,
        "plot(d[,1], type='l', lwd=4, ylim=c(min(d)-.5, max(d)+.5), main='%s')" % name,
        "lines(d[,2], col='red', lwd=4)",
        "legend('bottomleft', legend=c('promoter', 'enhancer'), col=c('black', 'red') , pch=20)",
        "abline(v=7, lty=2)",
        "abline(v=14, lty=2)",
        "abline(v=21, lty=2)",
        "dev.off",
        ]))
    f.close()
    os.system('/software/R/3.0.2/bin/Rscript %s' % (name + '.R'))
    os.system('rm %s %s' % (name + '.R', name))
    # os.system('convert %s.pdf %s.png' % (name, name))
    return (name + '.pdf')

    
def process(datafile, anchorfile, distalfile, headerfile):
    data = loadData(datafile)
    headerIndex = loadHeader(headerfile)
    S = calculateCorrelationMatrix(data)
    C = calculateCovarianceMatrix(data)
    S = C*S
    # C = 0.001*np.identity(len(MARKS))
    # N = calculateNormalization(data)
    M, SD = calculateMeanSdMatrix(datafile)
    offset = 2000000  ## search 4 Mbp up/downstream
    width_of_gaussian = 1000
    queryFile = pysam.Tabixfile(distalfile)
    count = 0
    with open(anchorfile) as inf:
        for record in csv.reader(inf, delimiter='\t'):
            X = convert2MatrixProximal(record, queryFile, width_of_gaussian, headerIndex, M, SD)
            X_slope = convertToSlope(X)
            start = int(record[ headerIndex['DHSpeak'] ]) - offset
            if start < 0:
                start = 0
            end = int(record[ headerIndex['DHSpeak'] ]) + offset
            detectedDistalRegions = [i for i in queryFile.fetch(record[ headerIndex['chrom'] ], start, end)]
            for detectedDistalRegion in detectedDistalRegions:
                Y = convert2MatrixDistal(detectedDistalRegion, queryFile, width_of_gaussian, headerIndex, M,  SD)
                Y_slope = convertToSlope(Y)
                count += 1
            # for xs in xrange(20):
            #     tmpX = []
            #     tmpY = []
            #     for i in xrange(len(MARKS)):
            #         tmpX.append(np.random.normal(0., 1., len(TIMES)))
            #         tmpY.append(np.random.normal(0., 1., len(TIMES)))

            #     X = np.array(tmpX)
            #     Y = np.array(tmpY)
            #     Z = np.copy(X)
            #     # adding the noise
            #     for i in xrange(len(MARKS)):
            #         X[i] += np.random.normal(0., C[i][i], len(TIMES))
            #         Y[i] += np.random.normal(0., C[i][i], len(TIMES))
            #         Z[i] += np.random.normal(0., C[i][i], len(TIMES))
                score = dependencyScore(X, Y, C, S)
                score_slope = dependencyScore(X_slope, Y_slope, C, S)
                score_flat = flatLikelihoodScore(X, Y, C, S)
                if score > 0 or True:
                    distalRecord = detectedDistalRegion.split()
                    if distalRecord[headerIndex['index']] == '':
                        distalRecord[headerIndex['index']] == 'unknown' 
                    print '\t'.join([
                        record[headerIndex['index']],
                        distalRecord[headerIndex['index']],
                        str( score ),
                        str( score_flat ),
                        str(score_slope),
                        str( int(record[1]) - int(distalRecord[1]) ),
                        record[headerIndex['Genename']],
                        ])
                if (score_flat>0.):
                    fname = record[headerIndex['index']] + '_' + distalRecord[headerIndex['index']]
                    f = open(fname, 'w')
                    for mark, m in zip(MARKS, range(len(MARKS))):
                        for t in xrange(len(TIMES)):
                            f.write('%f\t%f\n' % (X[m][t], Y[m][t]))
                    f.close()
                    plot(X, Y, fname)                                    

                
                

                        
def main():
    args = arguments()
    process(args.datafile, args.anchorfile, args.distalfile, args.headerfile)
    

if __name__ == '__main__':
    main()

    
