#! /software/bin/python

import sys
import os
import argparse
import re 

def arguments():
    parser = argparse.ArgumentParser(description="""
    given a directory containing BAM files, and also a directory of BED files that contains
    the interesting intervals, it submits jobs to Vital-IT for calculating the
    coverage of the BAM file over the genomic intervals. 
    """)
    parser.add_argument('-d', '--BAM', dest='bamDir', action='store',
                       type=str, required=True,
                       help='A directory that contains the BAM files.')
    parser.add_argument('-b', '--BED', dest='bedDir', action='store',
                       type=str, required=True,
                       help='BED files directory.') 
    parser.add_argument('-o', '--output', dest='destDir', action='store',
                       type=str, required=True,
                       help='The directory to store the outputs.')       
    args = parser.parse_args()
    return args


def submitJobs(bamFile, bedFile, dest):
    fname = os.path.basename(bamFile)
    bedFilename = os.path.basename(bedFile)
    bashFile = '_'.join([
        os.path.basename(bamFile),
        os.path.basename(bedFile) + '.sh',
        ])
    cmd = ' '.join([
        "bedtools",
        "coverage",
        '-abam "%s"'  % bamFile,
        '-b "%s"' % bedFile,
        ">",
        '"%s_%s"' % (os.path.join(dest, re.sub("\.bam$", "", fname)), bedFilename),
        ])
    with open(bashFile, 'w') as outf:
        outf.write('\n'.join( [
            "#!/bin/bash",
            "#BSUB -L /bin/bash",
            '#BSUB -o "%s_%s.stdout"' % (fname, bedFilename),
            '#BSUB -e "%s_%s.stderr"' % (fname, bedFilename),
            '#BSUB -J "%s_%s"' % (fname, bedFilename),
            '#BSUB -R "rusage[mem=50000]"',
            "#BSUB normal",
            "",
            "cd %s" % dest,
            cmd,
            ""
            ] ))
    os.system('bsub < %s' % bashFile )
    return bashFile            


def mergeResultsJobs(bamFile, dest):
    seq = 'bcdefghijklmnopqrst'
    fname = os.path.basename(bamFile)
    bashFile = os.path.basename(bamFile) + '.sh'
    cmd = ' '.join([
        "cat %s_xaa" % re.sub("\.bam$", "", fname),
        '|',
        'bedtools sort -i -',
        '> ',
        re.sub("\.bam$", ".bed", fname) +'\n',            
        ])    
    for s in seq:
        cmd += ' '.join([
            "cat %s_xa%s" % (re.sub("\.bam$", "", fname), s),
            '|',
            'bedtools sort -i -',
            '>>',
            re.sub("\.bam$", ".bed", fname) +'\n',            
            ])
    with open(bashFile, 'w') as outf:
        outf.write('\n'.join( [
            "#!/bin/bash",
            "#BSUB -L /bin/bash",
            '#BSUB -o "%s.stdout"' % (fname),
            '#BSUB -e "%s.stderr"' % (fname),
            '#BSUB -J "%s"' % (fname),
            '#BSUB -R "rusage[mem=50000]"',
            "#BSUB normal",
            "",
            "cd %s" % dest,
            cmd,
            ""
            ] ))
    os.system('bsub < %s' % bashFile )
    return bashFile 


def main():
    args = arguments()
    bamFiles = [os.path.join(args.bamDir, f) \
                for f in os.listdir(args.bamDir) if (re.search('\.bam$', f) and (not re.search('ZTall', f)))]
    bedFiles = [os.path.join(args.bedDir, f) \
                for f in os.listdir(args.bedDir)]
    for bamFile in bamFiles:
        mergeResultsJobs(bamFile, args.destDir)
    # for bamFile in bamFiles:
    #     for bedFile in bedFiles:
    #         submitJobs(bamFile, bedFile, args.destDir)
        

if __name__ == '__main__':
    main()
