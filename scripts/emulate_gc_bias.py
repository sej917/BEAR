#!/usr/bin/python
import zlib
import bz2
import numpy
import multiprocessing
from multiprocessing.managers import BaseManager, DictProxy
from collections import defaultdict, namedtuple, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from functools import partial
import sys, csv, StringIO, random, decimal, argparse
import ctypes as ct


#Return the GC content (as a percentage) of a given string
def gc_count(seq):
        return int(100 * ((seq.upper().count('G') + seq.upper().count('C') + 0.0) / len(seq)))


#Map the positions of all sequences with a specified length for all possible GC content values
def mapbuild(length, genome, max_length):
        genome_len = len(genome)
	print "working on length " + str(length)
        arr = numpy.zeros( (max_length+1,101))
        table = defaultdict(list)
        for start in range(0,genome_len):
		if(start % 1000000 == 0):
			print "at position " + str(start)
                gc = gc_count(genome[start:start+length])
                l = table[ (length, gc) ]
                l.append(start)
                table[ (length, gc) ] = l
                arr[length,gc] += 1
        return table

		
if __name__ == "__main__":

	
	parser = argparse.ArgumentParser(description='Generate reads according to a GC profile derived from a set of user-supplied NGS reads')
	parser.add_argument('-i', metavar='<genome_fasta>', dest="ref", help="Multi-FASTA file containing genomic sequences from which reads will be sampled.")
	parser.add_argument('-r', metavar='<read_file>', dest="reads", help="Multi-FASTA file containing reads from which to derive GC and read-length parameters.")
	parser.add_argument('-o', metavar='<output_file>', dest="output", help="Name for output file containing simulated uniform-length reads")
	parser.add_argument('-t', metavar='<total_reads>', type=int, dest="total", help="The total number of reads to sample from all genomes")
	parser.add_argument('-max', metavar="<max_read>", type=int, dest="max_read", help="Longest read length to simulate")
	parser.add_argument('-min', metavar="<min_read>", type=int, dest="min_read", help="Shortest read length to simulate")
	parser.add_argument('-p', metavar="<num_process>", type=int, dest="num_processes", help="Number of processes")
	parser.add_argument('-a', metavar='<abundance_file>', dest="abund", help="Tab-delimited abundance file with an abundance value for each corre- sponding genome sequence in <reference fasta>")
	args = parser.parse_args()

	if (args.max_read is not None) and (args.min_read is not None):
		assert args.max_read >= args.min_read
	total_reads = 0.0
	max_read = 0
	lens = []
	pairs = []
	for i in range(0,2000):
		lens.append(0)

	f1 = open(args.reads)
	print "Loading fastq file.."
	
	for read in SeqIO.parse(f1, 'fastq') :
		total_reads += 1.0
		readlen = len(read.seq)
		if(readlen > max_read):
			max_read = readlen
		lens[readlen] += 1.0
		if(readlen <= args.max_read) and (readlen >= args.min_read):
			pairs.append( (len(read.seq),gc_count(read.seq)) )
	f1.close()
	print "FASTQ file loaded"

	cond_dist = numpy.ndarray((max_read+1,101))
	
	f3 = open(args.abund);
	div_file = csv.reader(f3, delimiter='\t')
	species=[]
	diversity=[]

	for row in div_file:
	        species.append(row[0][1:])
       		diversity.append(decimal.Decimal(row[1]))	

        arr = range(args.min_read,args.max_read+1)
        arr2 = []
        total = 0.0
        for i in arr:
                total += 1.0/i
                if(total >= 0.01):
                        arr2.append(i)
                        total = 0
        print "Number of unique values to find: " + str(len(arr2))
        print "Processing genome..."


	f2 = open(args.ref)
	f4 = open(args.output, 'w')
	result = defaultdict(list)
	for i in SeqIO.parse(f2, "fasta"):
		genome_num=0
        	while(not(species[genome_num] in i.description)) :
        	        genome_num+=1
	        if(species[genome_num] in i.description) :
	                coverage=max(1, int((decimal.Decimal(diversity[genome_num])*args.total)))
			if(args.num_processes > 1):
				pool = multiprocessing.Pool(args.num_processes)
				partial_mapbuild = partial(mapbuild, genome=i.seq, max_length=args.max_read)
				newtable = pool.map(partial_mapbuild, arr2)
				pool.close()
				pool.join()
				for d in newtable: result.update(d)
			else:
				for readlen in arr2:
					newtable = mapbuild(readlen, i.seq, args.max_read)
					result.update(newtable)
	                limit=len(i.seq)
	                j=0
			
			for r in range(args.min_read, args.max_read+1):
				for gc in range(0, 101):
					if not r in arr2:
						result[(r,gc)] = result[(r-1,gc)]
			candidates = []
			
			print "Generating reads for " + i.description	
	                while j < coverage :
				params = random.choice(pairs)
				candidates = result[(params[0],params[1])]
				cur_count = len(candidates)
				k=1
				cur_abund = decimal.Decimal(diversity[genome_num])
				threshold = 0.01
				if(cur_abund < 0.01):
					threshold = cur_abund
				while(cur_count < threshold * args.total):
					#candidates = []
					candidates[:] = []
					for m in range(-k, k+1):
						for l in range(-k, k+1):
							key = (params[0]+m,params[1]+l)
							if key in result:
								candidates += result[(params[0]+m,params[1]+l)]
					cur_count = len(candidates)
					k+=1
				start = random.choice(candidates)
				end = start + params[0]
				header = str('>read' + str(j) + ' ' + i.description + ' ///' + ' positions ' + str(start) + ':' + str(end) + '\n')
				sequence = str(i.seq[start:end] + '\n')
				f4.write(header)
				f4.write(sequence)
				j+=1
			result.clear()
			for d in newtable: d.clear()
	f2.close()
	f4.close()

