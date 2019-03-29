#!/usr/bin/env bash
#The script for running DeepMod for single read.

python=/Users/andrewbailey/anaconda3/envs/embed_fast5/bin/python
DeepMod=/home/ubuntu/DeepMod/bin/DeepMod.py
workDir=/home/ubuntu/
readID=
refGenome=/home/ubuntu/refGenome/
base=
model=/home/ubuntu/DeepMod/train_mod/
PATH=$PATH:/home/ubuntu/minimap2/

${python} ${DeepMod} detect --wrkBase ${workDir} --Ref ${refGenome} --outFolder ${workDir} --Base ${base} --modfile ${model} --FileID ${readID} --alignStr minimap2 --threads 1
