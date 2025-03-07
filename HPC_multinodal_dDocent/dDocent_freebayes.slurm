#!/bin/bash

#SBATCH -J dDocent_freebayes                 # Name of the job
#SBATCH -o All_dDocent_freebayes.out       # Name of file that will have program output
#SBATCH -e All_dDocent_freebayes.err       # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 20                                # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal                               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=96:00:00

module load ddocent/2.7.8

ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | xargs -I {} -P 20 samtools view -@1 -b -h -1 -L mapped.{}.bed -o split.{}.bam cat-RRG.bam
ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | xargs -I {} -P 20 samtools index split.{}.bam
ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | xargs -I {} -P 20 freebayes -b split.{}.bam -t mapped.{}.bed -v raw.{}.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10 2>> fb.error.log
