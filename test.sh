# /bin/bash
#BSUB -L /bin/bash
#BSUB -o test.stdout
#BSUB -e test.stderr
#BSUB -J test 
#BSUB -q normal
#BSUB -f "/home/somidi/archive/CyclicX/ref > /scratch/cluster/monthly/somidi/Enhancer.Project/ref"
#BSUB -f "/home/somidi/archive/CyclicX/Dbp_bin_quantification.gz > /scratch/cluster/monthly/somidi/Enhancer.Project/Dbp_bin_quantification.gz"
#BSUB -f "/home/somidi/archive/CyclicX/Dbp_bin_quantification.gz.tbi > /scratch/cluster/monthly/somidi/Enhancer.Project/Dbp_bin_quantification.gz.tbi"

cd /scratch/cluster/monthly/somidi/Enhancer.Project
/home/somidi/codes/CyclicX/parse_file.py ref Dbp_bin_quantification.gz