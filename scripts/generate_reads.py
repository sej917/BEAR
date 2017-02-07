#!/usr/bin/python


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, csv, StringIO, random, decimal, argparse



parser = argparse.ArgumentParser(description='Generate uniform-length single or paired-end metagenomic reads.')
parser.add_argument('-r', metavar='<reference_fasta>', dest="ref", help="Multi-FASTA file containing genomic sequences from which reads will be sampled.")
parser.add_argument('-a', metavar='<abundance_file>', dest="abund", help="Tab-delimited abundance file with an abundance value for each corre- sponding genome sequence in <reference fasta>")
parser.add_argument('-o', metavar='<output_file>', dest="output", help="Name for output file containing simulated uniform-length reads")
parser.add_argument('-t', metavar='<total_reads>', type=int, dest="total", help="The total number of reads to sample from all genomes")
parser.add_argument('-l', metavar='<longest_read>', type=int, dest="length", help="The length, in bp, of the longest possible read to simulate")
parser.add_argument('-i', metavar='<insert_mean_length>', type=int, dest="insert", default="0", help="Average length of insert for paired-end reads.")
parser.add_argument('-s', metavar='<insert_stddev>', type=int, dest="stddev", default="0", help="Standard deviation of insert length for paired-end reads" )
parser.add_argument('-d', '--direction', action='store_true', dest="direction", help="Use this switch to generate reads in both forward and reverse orientations" )
args = parser.parse_args()


#Reference metagenome database file (FASTA)
f1 = open(args.ref);

#abundance file (tab-delimited .txt)
f2 = open(args.abund);

total_reads = args.total

max_read_length = args.length

insert_avg = args.insert
insert_stddev = args.stddev

if(insert_avg):
	f4 = open(args.output + '.1.fasta', 'w')
	f5 = open(args.output + '.2.fasta', 'w')
else:
	f4 = open(args.output, 'w')

frags=[]

div_file = csv.reader(f2, delimiter='\t')
species=[]
diversity=[]

lengths=[]
freqs=[]

for row in div_file:
	species.append(row[0][1:])
	diversity.append(decimal.Decimal(row[1]))


for i in SeqIO.parse(f1, 'fasta') :
	genome_num=0
	while(not(species[genome_num] in i.description)) :
		genome_num+=1
	if(species[genome_num] in i.description) :
		coverage=max(1, int((decimal.Decimal(diversity[genome_num])*total_reads)))

		limit=len(i.seq)
		for j in range(0, coverage) :
                	rand = random.random()
                	rand_length = 0
                	numLen = len(lengths)-1

<<<<<<< HEAD
			if( (insert_avg != 0) & (insert_stdev != 0)):
=======
			if( (insert_avg != 0) & (insert_stddev != 0)):
>>>>>>> 5cdada2dc2b92340e2e668b6273783525ebb29ad
				cur_insert = int(random.gauss(insert_avg, insert_stddev))
				if(limit > (max_read_length * 2 + cur_insert)):
					start1 = random.randint(0, limit-(2*max_read_length + cur_insert))
					end1 = start1 + max_read_length
					start2 = end1 + cur_insert
					end2 = start2 + max_read_length
				else:
					start1 = 0
					end1 = limit
					start2 = 0
					end2 = limit
				comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
				read1 = i.seq[start1:end1]
				read2 = ''.join([comp[b] for b in i.seq[end2:start2:-1]])
				if(args.direction):
					check = random.random()
					if(check < 0.5): #forward orientation
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read1)
						f5.write(">%s\n" % i.description)
                                		f5.write("%s\n" % read2)
					else: #reverse orientation
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read1[::-1])
						f5.write(">%s\n" % i.description)
						f5.write("%s\n" % read2[::-1])
			else:
				if(limit > max_read_length) :
					start=random.randint(0, limit-max_read_length)
					end=start+max_read_length
				else:
					start=0
					end=limit
				read = i.seq[start:end]
				if(args.direction):
					check = random.random()
					if(check < 0.5): #forward orientation
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read)
					else:
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read[::-1])

	if (genome_num >= len(species) ) :
		break;

f1.close()
f2.close()
f4.close()
