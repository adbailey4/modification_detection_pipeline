#The script for running minimap2 for single read.

minimap2=/home/ubuntu/minimap2/minimap2
samtools=/home/ubuntu/samtools/bin/samtools
refGenome=/home/ubuntu/refGenome/
readID=

${minimap2} -ax map-ont ${refGenome} ${readID}.fastq > ${readID}.sam
${samtools} view -Sb ${readID}.sam > ${readID}.bam
rm ${readID}.sam
