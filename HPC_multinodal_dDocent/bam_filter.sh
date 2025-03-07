#!/bin/bash

i=$1
echo "Filtering $i"
if [ ! -a $i.bam ] ; then bamtools filter -script $WORK/bin/filter.txt -in ${i}.bam -out ${i}-RG.bam
elif [ -a ${i}-RG.bam ] ; then bamtools filter -script $WORK/bin/filter.txt -in ${i}-RG.bam -out ${i}-tmp.bam; mv ${i}-tmp.bam ${i}-RG.bam
else echo "bam file for $i not found"
fi
