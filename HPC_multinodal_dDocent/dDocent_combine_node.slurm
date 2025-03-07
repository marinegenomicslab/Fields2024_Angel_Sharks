#!/bin/bash

#SBATCH -J dDocent_combine                      # Name of the job
#SBATCH -o dDocent_combine.out                  # Name of file that will have program output
#SBATCH -e dDocent_combine.err                  # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 20                                # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal,gpu                               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=24:00:00

module load ddocent/2.9.1
vcfcombine raw.*.vcf | sed -e 's/       \.\:/   \.\/\.\:/g' > cat.vcf
