#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o DBP_Longrange.stdout
#BSUB -e DBP_Longrange.stderr
#BSUB -J DBP
#BSUB -q normal

cd /scratch/cluster/monthly/somidi/DBP.Long.Range.Interaction
/home/somidi/codes/CyclicX/covarying_signal_identification.py -i Dbp -d Peaks_quantif_normalized_annotated.txt.distal.gz -g Peaks_quantif_normalized_annotated.txt.header > result