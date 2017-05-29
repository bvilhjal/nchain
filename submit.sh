#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 4
source /com/extra/BLAST/2.6.0/load.sh

python LD_pipeline.py cluster 'SM3_S10' &
python LD_pipeline.py cluster 'SM152B_S15' & 
python LD_pipeline.py cluster 'SM3_S10' &
wait