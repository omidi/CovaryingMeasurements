#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o Longrange.stdout
#BSUB -e Longrange.stderr
#BSUB -J LongRange
#BSUB -q normal

cd /scratch/cluster/monthly/somidi/Long.Range.Interaction
/home/somidi/codes/CyclicX/covarying_signal_identification.py -i Peaks_quantif_normalized_annotated.txt -d Peaks_quantif_normalized_annotated.txt.distal.gz -g Peaks_quantif_normalized_annotated.txt.header > result