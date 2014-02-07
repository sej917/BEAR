BEAR: Better Emulation for Artificial Reads
------------------------------------------
Created by Stephen Johnson, Brett Trost, Dr. Jeffrey R. Long, Dr. Anthony Kusalik
University of Saskatchewan, Department of Computer Science

BEAR is intended to be an easy-to-use collection of scripts for generating simulated WGS metagenomic reads with read lengths, quality scores, error profiles, and species abundances derived from real user-supplied WGS data.

BEAR is free for academic, non-commercial purposes. If you would like to use it for commercial purposes, please contact us.

BEAR is implemented as a collection of Perl and Python scripts, and is known to work with Perl v5.14.2 and Python v2.7.3. BEAR has the following dependencies:

    -Perl Getopt::Long, Bio::SeqIO, List::Util modules
    -Python Bio and Numpy packages
    -Python sys, csv, StringIO, random, decimal, argparse modules
    -DRISEE, which can be downloaded at https://github.com/MG-RAST/DRISEE

