#!/usr/bin/python
#Calculates all of the lengths or sequence from a multifasta file

from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
 print(output_line)


#from https://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/
#by GummyBear
