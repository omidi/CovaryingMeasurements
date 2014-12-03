#! /software/bin/python

import csv
import sys
import os
import numpy as np



def mleFit(rep1, rep2):
    N, M = np.sum(rep1), np.sum(rep2)
    logFreqDiffSqr = np.power(np.log(rep1/N) - np.log(rep2/M), 2.)
    countsInverseSum = 1./rep1 + 1./rep2
    sigmaLowerBound, sigmaUpperBound = 1.e-12, 1.e+3
    def calculateLikelihoodDerivative(sigma):
        sigmaSqr = 2.*sigma**2.
        likelihood = (logFreqDiffSqr / \
                        (2.*np.power((sigmaSqr + countsInverseSum), 2.)))
        likelihood -= 1.0 / (2*(sigmaSqr + countsInverseSum))
        return np.sum(likelihood)
    epsilon = 1.e-12
    while True:
        likelihoodLower = calculateLikelihoodDerivative(sigmaLowerBound)
        likelihoodUpper = calculateLikelihoodDerivative(sigmaUpperBound)
        sigmaMidd = (sigmaUpperBound + sigmaLowerBound) / 2.0
        # print sigmaLowerBound, sigmaMidd, sigmaUpperBound, likelihoodLower, likelihoodUpper
            
        likelihoodMidd = calculateLikelihoodDerivative(sigmaMidd)
        if likelihoodMidd > 0 and likelihoodLower > 0:
            sigmaLowerBound = sigmaMidd
        elif likelihoodMidd > 0 and likelihoodLower < 0:            
            sigmaUpperBound = sigmaMidd
        elif likelihoodMidd < 0 and likelihoodLower > 0:
            sigmaUpperBound = sigmaMidd
        else:
            sigmaLowerBound = sigmaMidd        
        # finishing condition
        if (sigmaUpperBound - sigmaLowerBound) < epsilon:
            break
    return sigmaMidd


def empiricalFit(rep1, rep2):
    N, M = np.sum(rep1), np.sum(rep2)
    logRatio = np.log(rep1/N) - np.log(rep2/M)
    meanRatio = np.mean(logRatio)
    sigmaEstimate = 1./(len(rep1)-1.) * np.sum(np.power(logRatio - meanRatio, 2.0))
    return sigmaEstimate
            

def main():
    if len(sys.argv)!=3:
        print "Usage: prog.py data_file mark"
        exit(1)
    
    mark = sys.argv[2]
    rep1_mark = mark + 'ZT02'
    rep2_mark = mark + 'ZT26'
    rep1, rep2 = [], []
    with open(sys.argv[1]) as inf:
        for record in csv.DictReader(inf, delimiter='\t'):
            if float(record[rep1_mark]) > 5. and float(record[rep2_mark]) >5.:
                    rep1.append(float(record[rep1_mark]))
                    rep2.append(float(record[rep2_mark]))
    sigmaEMP = empiricalFit(np.array(rep1), np.array(rep2))
    sigmaMLE = mleFit(np.array(rep1), np.array(rep2))
    print sigmaMLE, sigmaEMP    


if __name__ == '__main__':
    main()
