#! /software/bin/python

import sys
import os
import argparse
import re 

def arguments():
    parser = argparse.ArgumentParser(description="""
    given a directory containing BAM files, and also a BED file that contains
    the interesting intervals, it submits jobs to Vital-IT for calculating the
    coverage of the BAM file over the genomic intervals. 
    """)
    parser.add_argument('-d', '--BAM', dest='bamDir', action='store',
                       type=str, required=True,
                       help='A directory that contains the BAM files.')
    parser.add_argument('-b', '--BED', dest='bedFile', action='store',
                       type=str, required=True,
                       help='Input BED file that has genomic intervals.') 
    parser.add_argument('-o', '--output', dest='destDir', action='store',
                       type=str, required=True,
                       help='The directory to store the outputs.')       
    args = parser.parse_args()
    return args


def submitJobs(bamFile, bedFile, dest):
    fname = os.path.basename(bamFile)
    cmd = ' '.join([
        "bedtools",
        "coverage",
        '-abam "%s"'  % bamFile,
        '-b "%s"' % bedFile,
        ">",
        '"%s"' % os.path.join(dest, re.sub("\.bam$", "", fname)),
        ])
    with open(fname + '.sh', 'w') as outf:
        outf.write('\n'.join( [
            "#!/bin/bash",
            "#BSUB -L /bin/bash",
            '#BSUB -o "%s.stdout"' % fname,
            '#BSUB -e "%s.stderr"' % fname,
            '#BSUB -J "%s"' % fname,
            '#BSUB -R "rusage[mem=50000]"',
            "#BSUB normal",
            "",
            "cd %s" % dest,
            cmd,
            ""
            ] ))
    os.system('bsub < %s' % (fname + '.sh') )
    return fname + '.sh'            


def main():
    args = arguments()
    bamFiles = [os.path.join(args.bamDir, f) \
                for f in os.listdir(args.bamDir) if (re.search('\.bam$', f) and (not re.search('ZTall', f)))]
    for bamFile in bamFiles:
        submitJobs(bamFile, args.bedFile, args.destDir)
        

if __name__ == '__main__':
    main()
