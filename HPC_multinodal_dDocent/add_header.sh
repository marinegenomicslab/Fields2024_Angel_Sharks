#!/bin/bash

i=$1
echo "Checking $i"
if [ $(samtools view -H $i-RG.bam | head -n1 | grep HD | wc -l) -eq 0 ]; then 
echo "Adding HD heading to $i"
samtools view -h $i-RG.bam | sed -e '1 s/^/@HD\tVN:1.6\tSO:coordinate\n/' | samtools view -b > ${i}.bam; fi
