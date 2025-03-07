#!/bin/bash

#SBATCH -J dDocent_cat                       # Name of the job
#SBATCH -o dDocent_cat.out                 # Name of file that will have program output
#SBATCH -e dDocent_cat.err                 # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 20                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal                               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=96:00:00

module load ddocent/2.7.8

mkdir tmp
cat *.cov.stats | sort -k1,1 -k2,2n --temporary-directory=./tmp | bedtools merge -i - -c 4 -o sum > cov.stats

ls *-RG.bam > bamlist.list
samtools merge -@20 -b bamlist.list -f cat-RRG.bam
samtools index cat-RRG.bam
wait
bedtools merge -i cat-RRG.bam -bed >  mapped.bed
rm -r tmp
