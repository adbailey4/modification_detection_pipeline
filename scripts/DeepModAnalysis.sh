#!/usr/bin/env bash
#The script for merging DeepMod results.

python=/home/ubuntu/anaconda3/bin/python3.6
DeepMod=/home/ubuntu/DeepMod/tools/
workDir=/home/ubuntu/
refGenome=/home/ubuntu/refGenome/
base=C
runID=
chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y

${python} ${DeepMod}/sum_chr_mod.py ${workDir} ${base} ${runID} ${chromosomes}
#merge DeepMod results

${python} ${DeepMod}/generate_motif_pos.py ${refGenome} ${workDir} C CG 0
#5mC in CpG islands

${python} ${DeepMod}/hm_cluster_predict.py ${runID} ${workDir}
#detecting modified CpG island
