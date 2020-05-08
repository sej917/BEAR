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

if(args.direction):
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
comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 
    'N':'N', 'R':'Y', 'Y':'R', 'K':'M', 
    'M':'K', 'W':'W', 'S':'S', 'D':'H', 
    'H':'D', 'B':'V', 'V':'B', 'a':'t', 
    't':'a', 'g':'c', 'c':'g', 'n':'n', 
    'r':'y', 'y':'r', 'k':'m', 'm':'k', 
    'w':'w', 's':'s', 'd':'h', 'h':'d', 
    'b':'v', 'v':'b'}

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

			if(args.direction): ##Paired end
				cur_insert = int(random.gauss(insert_avg, insert_stddev))
				if(limit > cur_insert):
					start1 = random.randint(0, limit-cur_insert)
					end1 = start1 + max_read_length
					start2 = start1 + cur_insert - max_read_length
					end2 = start2 + max_read_length
				else:
					start1 = 0
					end1 = limit
					start2 = 0
					end2 = limit
				read1 = i.seq[start1:end1]
				read2 = ''.join([comp[b] for b in i.seq[end2:start2:-1]])
				f4.write(">%s\n" % i.description)
				f4.write("%s\n" % read1)
				f5.write(">%s\n" % i.description)
				f5.write("%s\n" % read2)
			else: 
				if(limit > max_read_length) :
					start=random.randint(0, limit-max_read_length)
					end=start+max_read_length
				else:
					start=0
					end=limit
				read = i.seq[start:end]
				if(args.direction and random.random() < 0.5):
                    #reverse orientation
					read = ''.join([comp[b] for b in i.seq[end:start:-1]])
					f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read)
					
				else:
                    #forward orientation
					f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read)

	if (genome_num >= len(species) ) :
		break;

f1.close()
f2.close()
f4.close()
