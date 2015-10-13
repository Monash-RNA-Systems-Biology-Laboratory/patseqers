#!/usr/bin/env python
from __future__ import print_function
import sys


"""Tool for clipping all copies of a single base off the start of FASTQ reads. By Andrew Pattison
Usage: python clip_first_bases.py base_to_clip Input_filename Output_filename"""


def clip_bases (base_to_clip, fqfile, output_filename):
	with open(output_filename, 'w') as f:
		number_off_clipped_bases = 0
		reads_clipped = 0
		total_reads = 0
		total_reads = len (fqfile) /4 
		count = 0
		for line in range(0,len(fqfile)):
			if line % 4 == 1 and fqfile[line][0] == base_to_clip:	
				reads_clipped = reads_clipped +1
				count = 0		
				for base in fqfile[line]:
					if base == base_to_clip:
						count = count + 1
					else:
						break
				number_off_clipped_bases = number_off_clipped_bases+count
				new_line = fqfile[line][count:]
				f.write(new_line)
				#print("clipping read",(line / 4)+1,"\n", count, "bases clipped")
				


			elif line % 4 == 3 and count != 0:
				new_qual_line = fqfile[line][count:]
				f.write(new_qual_line)
				count = 0
			else:
				f.write(fqfile[line])
	print("total clipped bases =", number_off_clipped_bases)
	print("total reads clipped =", reads_clipped)
	print("Percentage of reads clipped",  (float(reads_clipped) /float(total_reads) *100))



			




if __name__ == "__main__":

	base = sys.argv[1]

	fastq_file = sys.argv[2]

	output_filename = sys.argv[3]

	fqfile = open (fastq_file, "rU").readlines()

	clip_bases(base, fqfile, output_filename)