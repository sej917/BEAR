#!/usr/bin/env python

import os, sys, re, time, datetime, hashlib, shutil
import subprocess
import cStringIO
import string, random
from Bio import SeqIO
from random import randrange
from optparse import OptionParser
from multiprocessing import Pool

__doc__ = """
Script to calculate sequence error.
Input: fasta/fastq file to get error on
Output: error matrix file
STDOUT: Runtime summary stats"""

LOG_FILE = ''
ITER_MAX = 0
CONV_MIN = 0
PREF_LEN = 0

def ip (a,b, d):   # takes two sequences and offset and returns (number of matches, overlap length)
   assert len(b) >= len(a)
   if d < 0:
      r = len(a) + d
   elif (len(a) + d) > len(b):
      r =  len(b) - d
   else:
      r= len(a)
   c=0
   for i in range( 0, r ) :
      if (a[i] == b[i+d] and a[i] != "N") :
         c+=1
   if r==0:  r=1
   return c ,r

def align ( sequence, adapter ):  # takes two sequences, returns best alignment (number of matches, alignlength, offset)
   if len(sequence) > len(adapter) :
      a, b = adapter, sequence
   else:
      a, b = sequence, adapter
   assert len(b) >= len(a)
   la = len(a)
   lb = len(b)
   bestm = 0
   bestr = 0
   besto = 0
   for i in range(-len(a), len(b)):
      (m,r) = ip(a, b, i )
      c = float(m)/float(r)
      if m > bestm and c >= MINALIGNID and r >= MINOVERLAP :
         bestm = m
         besto = i
         bestr = r
   return (bestm , bestr, besto)

def bestalign( sequence, adapters):
   type(sequence)
   type(adapters)
   besta = 0
   bestk = 0
   for key in adapters.keys():
      (m,r,o) = align(sequence, adapters[key])
      if m > besta:
         besta = m
         bestk = key
   return (besta, bestk)

def write_file(text, fname, append=None):
    if append:
        mode = 'a'
    else:
        mode = 'w'
    outhdl = open(fname, mode)
    outhdl.write(text)
    outhdl.close()

def run_cmd(cmd, in_pipe=None, output=None):
    if not output:
        output = subprocess.PIPE
    if in_pipe:
        proc_in = subprocess.Popen( in_pipe, stdout=subprocess.PIPE )
        proc = subprocess.Popen( cmd, stdin=proc_in.stdout, stdout=output, stderr=subprocess.PIPE )
    else:
        proc = subprocess.Popen( cmd, stdout=output, stderr=subprocess.PIPE )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise IOError("%s\n%s"%(" ".join(cmd), stderr))
    return stdout, stderr

def random_truncate(items, cutoff):
    if (cutoff < 1) or (len(items) < 2) or (len(items) <= cutoff):
        return items
    # randomize array (fisher yates shuffle)
    i = len(items)
    while i > 1:
        i = i - 1
        j = randrange(i) # 0 <= j <= i-1
        items[j], items[i] = items[i], items[j]
    return items[:cutoff]

def random_str(size=6):
    chars = string.ascii_letters + string.digits
    return ''.join(random.choice(chars) for x in range(size))

# kk 11-9-12 # def seq_stats(in_file, fformat, verb):
#def seq_stats(in_file, opts.seq_type, verb):
def seq_stats(in_file, fformat, verb):
    #fout, ferr = run_cmd(['seq_length_stats.py', '-f', '-i', in_file, '-t', fformat])
    fout, ferr = run_cmd(['seq_length_stats.py', '-f', '-i', in_file, '-t', fformat])
    if LOG_FILE and (fout or ferr):
        write_file("\n".join([fout, ferr]), LOG_FILE, 1)
    if verb and (fout or ferr):
        sys.stdout.write("\n".join([fout, ferr]))
    lines = fout.strip().split('\n')
    stats = {}
    for l in lines:
        k, v = l.split('\t')
        stats[k] = v
    return stats

# kk 11-9-12 # def filter_seqs(in_file, out_file, stats, seqper, ambig_max, stdev_multi, filter_min, fformat):

#def filter_seqs(in_file, out_file, stats, seqper, ambig_max, stdev_multi, filter_min, opts.seq_type):
def filter_seqs(in_file, out_file, stats, seqper, ambig_max, stdev_multi, filter_min, fformat):
    # get stats
    avg_len = float(stats['average_length'])
    sdv_len = float(stats['standard_deviation_length'])
    min_len = int(avg_len - (sdv_len * stdev_multi))
    max_len = int(avg_len + (sdv_len * stdev_multi))
    if min_len < filter_min: min_len = filter_min + 1
    if max_len < min_len: max_len = min_len + 1
    # get filtered fasta file
    output_hdl = open(out_file, 'w')
    input_hdl = open(in_file, 'rU')
    new_num = 0
    try:
        # kk 11-9-12 # for rec in SeqIO.parse(input_hdl, fformat):
        for rec in SeqIO.parse(input_hdl, fformat):
            rnd_num = random.random()
            if seqper < rnd_num:
                continue
            ambig = 0
            slen = len(rec.seq)
            if (slen < min_len) or (slen > max_len):
                continue
            for char in rec.seq:
                if (char == 'n') or (char == 'N'): ambig += 1
            if ambig > ambig_max:
                continue
            new_num += 1
            output_hdl.write(">%s\n%s\n" % (rec.id, rec.seq)) #output_hdl.write(">%s\n%s\n" %(random_str(), rec.seq))
    finally:
        input_hdl.close()
        output_hdl.close()
    return new_num

def get_contaminated_md5_prefixes_from_fasta(in_file, prefix_len, adapters):
    contaminated_md5s = {}
    input_hdl = open(in_file, 'rU')
    try:
        for rec in SeqIO.parse(input_hdl, 'fasta'):
            seq = str(rec.seq).upper()
            prefix = seq[:prefix_len]
            md5 = hashlib.md5( prefix ).hexdigest()
            (a,b) = bestalign( rec.seq, adapters )
            if a > 0:
                contaminated_md5s[md5] = 1
    finally:
        input_hdl.close()
    return contaminated_md5s

def bin_replicate_seqs(in_file, out_file, tmp_dir, prefix_len, nodes):
    tmp_file = os.path.join(tmp_dir, os.path.basename(out_file)+'.tmp')
    tmp_hdl = open(tmp_file, 'w')
    input_hdl = open(in_file, 'rU')
    try:
        for rec in SeqIO.parse(input_hdl, 'fasta'):
            seq = str(rec.seq).upper()
            md5 = hashlib.md5( seq[:prefix_len] ).hexdigest()
            tmp_hdl.write("%s\t%s\n" %(md5, rec.id))
    finally:
        input_hdl.close()
        tmp_hdl.close()
    smem = str(nodes * 2 * 1024) + 'M'
    run_cmd(['sort', '-T', tmp_dir, '-S', smem, '-t', "\t", '-k', '1,1', '-o', out_file, tmp_file])
    #os.remove(tmp_file)

def filter_bins(fname, bin_min, total_max):
    bins = {}
    bin_num_out = open(fname+'.sum', 'w')
    run_cmd(['uniq', '-c'], ['cut', '-f1', fname], bin_num_out)
    bin_num_out.close()
    bin_num_in = open(fname+'.sum', 'rU')
    try:
        for line in bin_num_in:
            n, b = line.strip().split()
            if int(n) >= bin_min:
                bins[b] = n
    finally:
        bin_num_in.close()
    return dict([(x,'') for x in random_truncate(bins.keys(), total_max)])

def get_sub_fasta(ids, index_seq, seq_file, sub_fasta):
    # get sub fasta
    id_echo = ['echo']
    id_echo.extend(ids)
    stdo,_x = run_cmd(['cdbyank', index_seq, '-d', seq_file], id_echo)
    # get min length
    seq_lens = []
    str_hdl = cStringIO.StringIO(stdo)
    for rec in SeqIO.parse(str_hdl, 'fasta'):
        seq_lens.append( len(rec.seq) )
    str_hdl.close()
    # truncate all to min
    min_seq = min(seq_lens)
    out_hdl = open(sub_fasta, 'w')
    str_hdl = cStringIO.StringIO(stdo)
    for rec in SeqIO.parse(str_hdl, 'fasta'):
        out_hdl.write(">%s\n%s\n"%(rec.id, rec.seq[:min_seq]))
    str_hdl.close()
    out_hdl.close()

def process_bin(bin_id):
    bin_path = os.path.join(TMP_DIR, bin_id)
    os.mkdir(bin_path)
    cmd = ['run_find_steiner.pl','-i',bin_path+'.fasta','-o',bin_path+'.score','-l',bin_path+'.log','-t',bin_path,'--max_iter',str(ITER_MAX),'--min_conv',str(CONV_MIN), '-s', seq_name+'.uc']
    sto, ste = run_cmd(cmd)
    if LOG_FILE and (sto or ste):
        write_file("\n".join([sto, ste]), LOG_FILE, 1)
    #os.remove(bin_path+'.fasta')
    #shutil.rmtree(bin_path)
    return bin_id

def process_data(data, match, error):
    head = data.pop(0)
    bps = head.split("\t")
    for i, d in enumerate(data):
        if not d: continue
        counts = filter(lambda x: x != '', d.split("\t"))
        counts = map(lambda x: int(x), counts)
        if len(counts) != 6: continue
        if len(match) <= i: match.insert( i, dict([(x,0) for x in bps]) )
        if len(error) <= i: error.insert( i, dict([(x,0) for x in bps]) )
        max_bp = max(counts)
        for j, c in enumerate(counts):
            if c == max_bp: match[i][ bps[j] ] += c
            elif c < max_bp: error[i][ bps[j] ] += c
    return bps, match, error

def create_output(bps, match, error, per):
    total = 0
    errs = dict([(x,0) for x in bps])
    stext = [ "#\t" + "\t".join(bps) + "\t" + "\t".join(bps) ]
    
    for i in range( len(match) ):
        rsum = 0
        row = []
        for b in bps:
            row.append( match[i][b] )
            rsum += match[i][b]
            if i > (PREF_LEN - 1):
                total += match[i][b]
        for b in bps:
            row.append( error[i][b] )
            rsum += error[i][b]
            if i > (PREF_LEN - 1):
                total += error[i][b]
                errs[b] += error[i][b]
        if per:
            row = map(lambda x: "%f"%(((x * 1.0) / rsum) * 100), row)
        else:
            row = map(lambda x: str(x), row)
        stext.append( "%d\t"%(i+1) + "\t".join(row) )

    err_head = map(lambda x: "%s_err"%x, bps)
    err_nums = map(lambda x: "%f"%(((errs[x] * 1.0) / total) * 100), bps)
    err_sum = ((sum(errs.values()) * 1.0) / total) * 100
    err_head.append('bp_err')
    err_nums.append("%f"%err_sum)

    stext.insert(0, "#\t"+"\t".join(err_nums))
    stext.insert(0, "#\t"+"\t".join(err_head))
    return err_sum, "\n".join(stext)


usage = "usage: %prog [options] input_seq_file output_stat_file\n" + __doc__
version = "%prog 1.5"

def main(args):
    global TMP_DIR, LOG_FILE, ITER_MAX, CONV_MIN, PREF_LEN, MINALIGNID, MINOVERLAP
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-p", "--processes", dest="processes", type="int", default=8, help="Number of processes to use [default '8']")
    #parser.add_option("-t", "--seq_type", dest="seq_type", default='fasta', help="Sequence type: fasta, fastq [default 'fasta']")
    parser.add_option("-t", "--seq_type", dest="seq_type", default='fasta', help="Sequence type: fasta, fastq [default 'fasta']")
    parser.add_option("-f", "--filter_seq", dest="filter", action="store_true", default=False, help="Run sequence filtering, length and ambig bp [default off]")
    parser.add_option("-r", "--replicate_file", dest="rep_file", default=None, help="File with sorted replicate bins (bin_id, seq_id) [default to calculate replicates]")
    parser.add_option("-d", "--tmp_dir", dest="tmpdir", default="/tmp", help="DIR for intermediate files (must be full path). Specified directory must already exist. Files are automatically deleted at analysis completion. [default '/tmp']")
    parser.add_option("-l", "--log_file", dest="logfile", default=None, help="Detailed processing related stats [default '/dev/null']")
    parser.add_option("", "--percent", dest="percent", action="store_true", default=False, help="Additional output (output_stat_file.per) with percent values [default off]")
    parser.add_option("", "--prefix_length", dest="prefix", type="int", default=50, help="Prefix length for replicate bins [default 50]")
    parser.add_option("-s", "--seq_max", dest="seq_max", type="int", default=1000000, help="Maximum number of reads to process (chosen randomly) [default 1000000]")
    parser.add_option("-a", "--ambig_bp_max", dest="ambig_max", type="int", default=0, help="Maximum number of ambiguity characters before rejection [default 0]")
    parser.add_option("-m", "--stdev_multiplier", dest="stdev_multi", type="float", default=2.0, help="Multiplier to stddev to get min and max seq lengths [default 2.0]")
    parser.add_option("-n", "--bin_read_min", dest="read_min", type="int", default=20, help="Minimum number of reads in bin to be considered [default 20]")
    parser.add_option("-x", "--bin_read_max", dest="read_max", type="int", default=1000, help="Maximum number of reads in bin to process (chosen randomly) [default 1000]")
    parser.add_option("-b", "--bin_num_max", dest="num_max", type="int", default=1000, help="Maximum number of bins to process (chosen randomly) [default 1000]")
    parser.add_option("-i", "--iter_max", dest="iter_max", type="int", default=10, help="Maximum number of iterations if alignment does not converge [default 10]")
    parser.add_option("-c", "--converge_min", dest="conv_min", type="int", default=3, help="Minimum number of iterations to identify convergence [default 3]")
    parser.add_option("-j", "--check_contam", dest="check_contam", action="store_true", default=False, help="Separate results for seqs with adapter contamination [default off]")
    parser.add_option("-o", "--minoverlap", dest="MINOVERLAP", type="int", default=10, help="Minimum overlap paramter for identifying adapter contamination [default 10]")
    parser.add_option("-e", "--fractionid", dest="MINALIGNID", type="float", default=0.9, help="Minimum alignment id for identifying adapter contamination [default 0.9]")
    parser.add_option("-g", "--database", dest="database", default="adapterDB.fna", help="Database fasta of adapter sequences [default adapterDB.fna]")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Write runtime summary stats to STDOUT [default off]")

    start_time = time.time()
    (opts, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("Incorrect number of arguments")
        
    # check inputs
    global seq_name
    (in_seq, out_stat) = args
    seq_name = in_seq
    if not (os.path.isfile(in_seq) and os.path.isdir(opts.tmpdir)):
       #parser.error("Invalid input files and/or tmp dir")
       error_string = ""
       if not (os.path.isfile(in_seq)): error_string = ( "Invalid input file (" + in_seq + ")\t;" )
       if not (os.path.isdir(opts.tmpdir)): error_string = ( error_string + "Invalid tmp dir (" + opts.tmpdir + ")\n" )
       parser.error(error_string)
    if opts.processes < 1: opts.processes = 1
    if opts.ambig_max < 0: opts.ambig_max = 0
    if opts.stdev_multi <=0: opts.stdev_multi = 2
    if opts.read_min > opts.read_max:
        parser.error("bin_read_min (%d) can not be greater than bin_read_max %d)"%(opts.read_min, opts.out_name))
    if opts.read_min < 1: opts.read_min = 1
    if opts.read_max < 2: opts.read_max = 2
    if opts.num_max < 1: opts.num_max = 1
    if opts.iter_max < 1: opts.iter_max = 1
    if opts.conv_min < 1: opts.conv_min = 1
    if opts.seq_max < 1: opts.seq_max = 2
    if opts.prefix < 10: opts.prefix = 10
    TMP_DIR = os.path.join(opts.tmpdir, random_str(8)+'.drisee')
    LOG_FILE = opts.logfile
    ITER_MAX = opts.iter_max
    CONV_MIN = opts.conv_min
    PREF_LEN = opts.prefix
    MINOVERLAP = opts.MINOVERLAP
    MINALIGNID = opts.MINALIGNID
    DBFILE = opts.database

    removeambigtable= string.maketrans('RYWSMKHBVDNrywsmkhbvdn','N'* 22)
    adapters = {}
    adapters[0] = ""
    for seq_record in SeqIO.parse(DBFILE, "fasta"):
        adapters[seq_record.description] =        string.upper(str(seq_record.seq).translate(removeambigtable))
        adapters["%s.R"%seq_record.description] = string.upper(str(seq_record.seq.reverse_complement()).translate(removeambigtable))

    os.mkdir(TMP_DIR)
    
    if opts.verbose: sys.stdout.write("Version:\t%s\n"%version)
    # seq stats
    stats = seq_stats(in_seq, opts.seq_type, opts.verbose)
    seqnum = int(stats['sequence_count'])
    seqper = float(opts.seq_max) / seqnum
    seqmax = 0
    if opts.verbose: sys.stdout.write("DRISEE will be run on %d of %d sequences\n"%(min(opts.seq_max, seqnum),seqnum))

    # random subselect / length filter
    if opts.filter:
        if opts.verbose: sys.stdout.write("Filtering with max ambig %d and stddev range x%f ... " %(opts.ambig_max,opts.stdev_multi))
        filter_file = os.path.join(TMP_DIR, os.path.basename(in_seq)+'.filter.fasta')
        seqmax = filter_seqs(in_seq, filter_file, stats, seqper, opts.ambig_max, opts.stdev_multi, opts.prefix, opts.seq_type)
        in_seq = filter_file
        opts.seq_type = 'fastq'
        if opts.verbose: sys.stdout.write("Done, %s sequences kept\n"%seqmax)
    # random subselect / uniquify seqs
    else:
        if opts.verbose: sys.stdout.write("Making sure seq ids are unique ... ")
        out_file = os.path.join(TMP_DIR, os.path.basename(in_seq)+'.uniq.fasta')
        out_hdl = open(out_file, 'w')
        in_hdl = open(in_seq, 'rU')
        try:
            for rec in SeqIO.parse(in_hdl, opts.seq_type):
                rnd_num = random.random()
                if seqper >= rnd_num:
                    seqmax += 1
                    out_hdl.write(">%s\n%s\n" %(rec.id, rec.seq)) #out_hdl.write(">%s\n%s\n" %(random_str(), rec.seq))
        finally:
            in_hdl.close()
            out_hdl.close()
        in_seq = out_file
        opts.seq_type = 'fastq'
        if opts.verbose: sys.stdout.write("Done, %s sequences kept\n"%seqmax)

    ### seq file is always fasta from here on
    # dereplication
    if not (opts.rep_file and os.path.isfile(opts.rep_file)):
        if opts.verbose: sys.stdout.write("Creating replicate bins, prefix size %d bps ... " %opts.prefix)
        rep_file = os.path.join(TMP_DIR, os.path.basename(in_seq)+'.derep')
        bin_replicate_seqs(in_seq, rep_file, TMP_DIR, opts.prefix, opts.processes)
        opts.rep_file = rep_file
        if opts.verbose: sys.stdout.write("Done\n")

    # index fasta file
    if opts.verbose: sys.stdout.write("Creating cdb index ... ")
    index_seq = os.path.join(TMP_DIR, os.path.basename(in_seq)+'.cidx')
    iout, ierr = run_cmd(['cdbfasta', in_seq, '-o', index_seq])
    if opts.logfile and (iout or ierr):
        write_file("\n".join([iout, ierr]), opts.logfile, 1)
    if opts.verbose: sys.stdout.write("Done\n")

    # filter bin set
    if opts.verbose: sys.stdout.write("Getting %d random bins with >= %d reads ... " %(opts.num_max,opts.read_min))
    bins = filter_bins(opts.rep_file, opts.read_min, opts.num_max)
    size = len(bins)
    if size == 0:
        msg = "No available bins >= %d to process, quiting\n"%opts.read_min
        if opts.logfile: write_file(msg, opts.logfile, 1)
        if opts.verbose: sys.stdout.write(msg)
        write_file('', out_stat)
        #shutil.rmtree(TMP_DIR)
        return 0
    if opts.verbose: sys.stdout.write("Done, %s bins found\n"%size)

    # check for adapter contamination in one sequence for each md5 in filtered set
    if opts.check_contam:
        # bins_to_rep_ids will contain a dictionary of each md5 in filtered set to one sequence of the sequence ids
        bins_to_rep_ids = {}
        # bins_to_seq_count will contain a dictionary of each md5 in filtered set to a count of the number of seqs in that bin
        bins_to_seq_count = {}
        dhdl = open(opts.rep_file, 'rU')
	try:
	    for line in dhdl:
		(bid, sid) = line.split()
                if bid in bins:
                    bins_to_rep_ids[bid] = sid
                    if bid in bins_to_seq_count:
                        bins_to_seq_count[bid] += 1
                    else:
                        bins_to_seq_count[bid] = 1
	finally:
	    dhdl.close()
        # generate sequence file for list of sequence ids
        rep_seqs_fasta = os.path.join(TMP_DIR, "rep_seqs.fasta")
        get_sub_fasta(bins_to_rep_ids.values(), index_seq, in_seq, rep_seqs_fasta)
        contaminated_md5s = get_contaminated_md5_prefixes_from_fasta(rep_seqs_fasta, opts.prefix, adapters)
        #os.remove(rep_seqs_fasta)

    # create trimmed bin fasta files
    to_process = []
    total_ids = 0
    ids = []
    curr = ''
    dhdl = open(opts.rep_file, 'rU')
    if opts.verbose: sys.stdout.write("Creating %d bin files with >= %d reads ..." %(size,opts.read_min))
    try:
        for line in dhdl:
            (bid, sid) = line.split()
            if not (bid and sid and (bid in bins)):
                continue
            if curr == '':
                curr = bid
            if bid != curr:
                if len(ids) > opts.read_max:
                    ids = random_truncate(ids, opts.read_max)
                bin_fasta = os.path.join(TMP_DIR, curr+".fasta")
                get_sub_fasta(ids, index_seq, in_seq, bin_fasta) ## seqs truncated to min
                if os.path.isfile(bin_fasta):
                    to_process.append(curr)
                    total_ids += len(ids)
                curr = bid
                ids = []
            ids.append(sid)
        if len(ids) > opts.read_max:
            ids = random_truncate(ids, opts.read_max)
        bin_fasta = os.path.join(TMP_DIR, curr+".fasta")
        get_sub_fasta(ids, index_seq, in_seq, bin_fasta) ## seqs truncated to min
        if os.path.isfile(bin_fasta):
            to_process.append(curr)
            total_ids += len(ids)
    finally:
        dhdl.close()
    if opts.verbose: sys.stdout.write("Done\n")

    # process bins
    min_proc = 1
    if opts.processes > min_proc:
        min_proc = opts.processes
    if len(to_process) < min_proc:
        min_proc = len(to_process)
    if opts.verbose: sys.stdout.write("Processing %d bins (%d sequences total) using %d threads ... "%(size,total_ids,min_proc))
    pool = Pool(processes=min_proc)
    finish = pool.map(process_bin, to_process,1)
    pool.close()
    pool.join()
    if opts.verbose: sys.stdout.write("Done\n")

    # merge results
    if len(finish) == 0:
        msg = "No bins were processed, quiting\n"
        if opts.logfile: write_file(msg, opts.logfile, 1)
        if opts.verbose: sys.stdout.write(msg)
        write_file('', out_stat)
        return 0
    bases = []
    match = []
    error = []
    # if check_contam is on, create stats for contaminated and non-contaminated seqs in addition to stats for all seqs
    if opts.check_contam:
        cbases = []
        cmatch = []
        cerror = []
        ncbases = []
        ncmatch = []
        ncerror = []
    if opts.verbose: sys.stdout.write("Merging scores from %d bins ... "%len(finish))
    for bid in finish:
        bin_score = os.path.join(TMP_DIR, bid+'.score')
        bin_log = os.path.join(TMP_DIR, bid+'.log')
        if opts.logfile and os.path.isfile(bin_log):
            lhdl = open(bin_log, 'rU')
            write_file(lhdl.read(), opts.logfile, 1)
            lhdl.close()
        if os.path.isfile(bin_score):
            shdl = open(bin_score, 'rU')
            data = shdl.read().split("\n")
            data_copy = list(data)
            bases, match, error = process_data(data_copy, match, error)
            if opts.check_contam:
                if bid in contaminated_md5s:
                    data_copy = list(data)
                    cbases, cmatch, cerror = process_data(data_copy, cmatch, cerror)
                else:
                    data_copy = list(data)
                    ncbases, ncmatch, ncerror = process_data(data_copy, ncmatch, ncerror)
            shdl.close()
    if opts.verbose: sys.stdout.write("Done\n")
    if opts.check_contam:
        contam = len(contaminated_md5s)
        ncontam = len(finish) - len(contaminated_md5s)
        contam_seq_count = 0
        ncontam_seq_count = 0
        for md5 in bins_to_seq_count:
            if md5 in contaminated_md5s:
                contam_seq_count += bins_to_seq_count[md5]
            elif md5 in finish:
                ncontam_seq_count += bins_to_seq_count[md5]

    err_score, score_text = create_output(bases, match, error, 0)
    write_file(score_text, out_stat)
    if opts.check_contam:
        if contam > 0:
            cerr_score, cscore_text = create_output(cbases, cmatch, cerror, 0)
            write_file(cscore_text, out_stat+".contaminated")
        if ncontam > 0:
            ncerr_score, ncscore_text = create_output(ncbases, ncmatch, ncerror, 0)
            write_file(ncscore_text, out_stat+".non_contaminated")

    if opts.percent:
        _err, per_text = create_output(bases, match, error, 1)
        write_file(per_text, out_stat+'.per')
        if opts.check_contam:
            if contam > 0:
                _err, cper_text = create_output(cbases, cmatch, cerror, 1)
                write_file(cper_text, out_stat+'.contaminated.per')
            if ncontam > 0:
                _err, ncper_text = create_output(ncbases, ncmatch, ncerror, 1)
                write_file(ncper_text, out_stat+'.non_contaminated.per')
    #print match
    #print bases
    #print error
    # cleanup
    #shutil.rmtree(TMP_DIR)
    end_time = time.time() - start_time
    if opts.verbose: sys.stdout.write("Completed in %s\n" %str(datetime.timedelta(seconds=end_time)))
    if opts.verbose: sys.stdout.write("Input seqs\t%d\nProcessed bins\t%d\nProcessed seqs\t%d\nDrisee score\t%f\n"%(seqmax,size,total_ids,err_score))
    if opts.verbose and opts.check_contam and contam > 0: sys.stdout.write("\nContam bins\t%d\nContam seqs\t%d\nDrisee score\t%f\n"%(contam,contam_seq_count,cerr_score))
    if opts.verbose and opts.check_contam and ncontam > 0: sys.stdout.write("\nNon-contam bins\t%d\nNon-contam seqs\t%d\nDrisee score\t%f\n"%(ncontam,ncontam_seq_count,ncerr_score))
    return 0
    

if __name__ == "__main__":
    sys.exit(main(sys.argv))
