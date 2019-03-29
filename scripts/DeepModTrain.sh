#The script for training DeepMod.

python=/home/ubuntu/anaconda3/bin/python3.6
DeepMod=/home/ubuntu/DeepMod/bin/DeepMod.py

${python} ${DeepMod} train --wrkBase umr --wrkBase2 sss --FileID mod_train --outFolder ./mod_output/train1
